#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.bed_files = "${baseDir}/input_bed"
params.ld_bfile = "$baseDir/input_bed/all_chr/all_genomes_09112019_all_chr"
params.phenotype = "${baseDir}/data/train_phenotype.tsv"
params.metabolites = "${baseDir}/data/met_indices.csv"
params.output = "${baseDir}/data"
params.grm_sparsity = 0.05
params.gwas_pval = 5.29e-11
params.filter_pval = 1e-5
params.max_threads = 12

def helpMessage() {
    log.info"""
    ~~~ GWAS workflow ~~~
    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf --resume
    A run with all parameters set would look like:
    > nextflow run main.nf --data_dir=./dat--refs=/my/references --single_end=false \\
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
    val(bed_files)

    output:
    path("gwas_chromosomes.txt")

    """
    echo "${bed_files.join('\n')}" > gwas_chromosomes.txt
    """
}

process build_grm {
    cpus params.max_threads

    input:
    path(files)
    path(gwas_chromosomes)

    output:
    path("sp_grm.*.*")

    """
    gcta64 --mbfile ${gwas_chromosomes} \\
        --make-grm --sparse-cutoff ${params.grm_sparsity} \\
        --thread-num ${task.cpus} --out sp_grm
    """

}

process gwas {
    publishDir "${params.output}/gwas_results/", mode: "copy", overwrite: false

    cpus 4

    input:
    tuple val(mid), val(metabolite)
    path(gwas_chromosomes)
    path(grm)
    path(genome_files)

    output:
    tuple val("${metabolite}"), path("gwas_${metabolite}.fastGWA.gz")

    """
    gcta64 --mbfile ${gwas_chromosomes} \\
        --grm-sparse sp_grm --fastGWA-mlm \\
        --pheno ${params.phenotype} --mpheno ${mid} \\
        --thread-num ${task.cpus} --out gwas_${metabolite} || touch gwas_${metabolite}.fastGWA
    gzip gwas_${metabolite}.fastGWA
    """
}

process filter_gwas {
    publishDir "${params.output}/gwas_filtered/", mode: "copy", overwrite: false

    cpus 1

    input:
    tuple val(metabolite), path(gwas)

    output:
    tuple val("${metabolite}"), path("filtered_${metabolite}.tsv.gz")

    """
    zcat ${gwas} \\
        | awk '{ if (\$10 < ${params.filter_pval} || \$10 == "P") { print }}' \\
        | gzip > filtered_${metabolite}.tsv.gz
    """
}

process remove_ld_first {
    cpus 2

    input:
    tuple val(metabolite), path(gwas_result)

    output:
    tuple val(metabolite), path("${metabolite}_pass1.clumped")

    """
    plink --bfile ${params.ld_bfile} \\
         --clump ${gwas_result} --clump-kb 250 \\
         --clump-p1 ${params.gwas_pval} --clump-p2 5e-8 \\
         --clump-r2 0.8 --maf 0.01 --hwe 1e-10 midp \\
         --out ${metabolite}_pass1 || \\
         echo "no associations, skipping"
    touch ${metabolite}_pass1.clumped
    """
}

process remove_ld_second {
    publishDir "${params.output}/clumped"

    cpus 2

    input:
    tuple val(metabolite), path(gwas_result), path(clumped)

    output:
    path("${metabolite}_pass2.clumped")

    """
    awk -F ' ' '{print \$3}' ${clumped} > snps.txt
    plink --bfile ${params.ld_bfile} \\
         --clump ${gwas_result} --extract snps.txt --clump-kb 250 \\
         --clump-p1 ${params.gwas_pval} --clump-p2 5e-8 \\
         --clump-r2 0.3 --maf 0.01 --hwe 1e-10 midp \\
         --out ${metabolite}_pass2 || \\
         echo "no associations, skipping"
    touch ${metabolite}_pass2.clumped
    """
}

process merge_results {
    publishDir "${params.output}", mode: "copy", overwrite: true

    cpus 4

    input:
    path(clumped)
    path(gwas_results)

    output:
    path("final_results.csv")

    script:
    """
    #!/usr/bin/env python

    from os import path
    import pandas as pd

    inputs = "${clumped}".split()

    results = []
    for res in inputs:
        metabolite = res.split("_pass")[0]
        try:
            r = pd.read_csv(res, sep="\s+", header=0)
            g = pd.read_csv(f"filtered_{metabolite}.tsv.gz", sep="\t", header=0)
        except pd.errors.EmptyDataError:
            continue
        print(f"{res} has significant SNPs")
        merged = r.merge(g[["SNP", "BETA", "SE", "A1", "A2"]], on="SNP")
        merged["metabolite"] = metabolite
        results.append(merged)
    results = pd.concat(results)
    results.to_csv("final_results.csv", index=False)
    """
}

workflow {
    Channel
        .fromPath("${params.bed_files}/*.bed")
        .ifEmpty { error "Cannot find any BED files in ${params.bed_files}!" }
        .map{row -> tuple(row.baseName.split("\\.bed")[0])}
        .set{bed}
    Channel
        .fromPath("${params.bed_files}/*.{bed,bim,fam}")
        .ifEmpty { error "Cannot find any files in ${params.bed_files}!" }
        .set{genome_files}
    Channel
        .fromPath("${params.metabolites}")
        .splitCsv(header: true)
        .set{mets}


    bed.collect() | list_data_by_chromosomes
    build_grm(genome_files.collect(), list_data_by_chromosomes.output)

    gwas(
        mets,
        list_data_by_chromosomes.output,
        build_grm.output,
        genome_files.collect()
    ) | filter_gwas

    remove_ld_first(filter_gwas.output)
    remove_ld_second(filter_gwas.output.join(remove_ld_first.output))

    merge_results(remove_ld_second.out.collect(), filter_gwas.output.map{row -> row[1]}.collect())
}
