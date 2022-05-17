#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.bed_files = "${baseDir}/input_bed"
params.ld_bfile = "$baseDir/input_bed/all_chr/all_genomes_09112019_all_chr"
params.phenotype = "${baseDir}/data/train.tsv"
params.metabolites = "${baseDir}/data/mets.csv"
params.output = "${baseDir}/data"
params.grm_sparsity = 0.05
params.gwas_pval = 5.29e-11

def helpMessage() {
    log.info"""
    ~~~ GWAS workflow ~~~
    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf --resume
    A run with all parameters set would look like:
    > nextflow run main.nf --data_dir=./data --single_end=false --refs=/my/references --single_end=false \\
                           --trim_front=5 --min_length=50 --quality_threshold=20 --read_length=150 --threshold=10
    General options:
      --bed_files                   Folder containing the BED files for the analysis.
      --ld_bfile                    Basepath of BED/BIM files used for LD corrections. Can be a joined
                                    version of the input bed files.
      --output [str]                The main data directory for the analysis (must contain `raw`).
      --threads [int]               The maximum number of threads a single process can use.

    GWAS options:
      --phenotype                   A tab-separated file containing the confounder-corrected phenotype data.
      --metabolites                 A text file containing a list of phenotype columns, one per line.
      --grm_sparsity                How sparse the GRM matrix should be.
      --gwas_pval                   P-value cutoff that denotes genome-wide significance.
    """.stripIndent()
}

params.help = false
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


process list_data_by_chromosomes {
    cpus 1

    input:
    path(bed_files)

    output:
    path("gwas_chromosomes.txt")

    """
    echo "${bed_files.join('\n')}" > gwas_chromosomes.txt
    """
}

process build_grm {
    cpus ${params.max_threads}

    input:
    path(gwas_chromosomes)

    output:
    path("sp_grm")

    """
    gcta64 --mbfile ${gwas_chromosomes} \\
        --make-grm --sparse-cutoff ${params.grm_sparsity} \\
        --thread-num ${task.cpus} --out sp_grm
    """

}

process gwas {
    publishDir "${params.output}/gwas_results/", mode: "copy", overwrite: True

    cpus 4

    input:
    tuple val(metabolite), path(gwas_chromosomes), path(grm)

    output:
    tuple val(metabolite_id), path("gwas_{metabolite_id}.tsv")

    """
    gcta64 --mbfile ${gwas_chromosomes} \\
        --grm-sparse ${grm} --fastgwa-mlm \\
        --pheno ${params.phenotype} --mpheno ${metabolite} \\
        --thread-num ${task.cpus} --out
    """
}

process remove_ld {
    publishDir "${params.output}/clumped"

    cpus 4

    input:
    tuple val(metabolite_id), path(gwas_result), path(gwas_chromosome), val(r2_cutoff)

    output:
    path("clumped_${metabolite_id}_${r2_cutoff}.tsv")

    """
    plink --bfile ${params.ld_bfile} \\
         --clump ${gwas_result} --clump-kb 250 \\
         --clump-p1 ${params.gwas_pval} --clump-p2 5e-8 \\
         --clump-r2 ${r2_cutoff} --maf 0.01 --hwe 1e-10 midp \\
         --out clumped_${metabolite_id}_${r2_cutoff}.tsv
    """
}

proces merge_results {
    publishDir "${params.output}", mode: "copy", overwrite: True

    cpus 4

    input:
    path(gwas_results)

    output:
    path("final_results.csv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    inputs = "${gwas_results}".split()

    results = []
    for res in inputs:
        r = pd.read_csv(res, sep="\t", header=True)
        results.append(r[r.P < ${params.gwas_pval}])
    results = pd.concat(results)
    results.to_csv("final_results.csv", index=False)
    """
}

workflow main {
    Channel
        .fromPath("${params.bed_files}/*.bed")
        .ifEmpty { error "Cannot find any BED files in ${params.bed_files}!" }
        .map{row -> tuple(row.baseName.split("\\.fastq")[0])}
        .set{bed}
    Channel
        .fromPath("${params.metabolites}")
        .splitText()
        .set{mets}


    bed | list_data_by_chromosomes | build_grm
    gwas_runs = mets.combine(list_data_by_chromosomes.output).combine(build_grm.output)
    gwas(gwas_runs)
    remove_ld(gwas.output.combine(list_data_by_chromosomes.output).combine(0.6))
    // remove_ld(gwas.output.combine(list_data_by_chromosomes.output).combine(0.1))

    merge_results(remove_ld.output)
}
