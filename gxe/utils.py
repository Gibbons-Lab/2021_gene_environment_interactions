"""Some general utility functions."""

from Bio import Entrez
import pandas as pd
import re
import seaborn as sns

Entrez.email = "cdiener@isbscience.org"


def rename_numerical_id(ID):
    """Rename a numerical ID to a valid column name.
    
    Parameters
    ----------
    ID : str
        The ID to adjust.
        
    Returns
    -------
    str
        The new ID.
    """
    if re.match("^[0-9]", str(ID)) != None:
        return 'X' + str(ID)
    else:
        return str(ID)

    
def rsid2gene(ids):
    """Find the genes for a set of rsids.
    
    Parameters
    ----------
    ids : array or list of str
        The SNV IDs to be looked up in dbSNP.
    
    Returns
    -------
    pandas.DataFrame
        A data frame containing the identified variants and 
        associated information such as genes, clinical significance,
        chromosome and functional class.
    """
    id_str = ",".join([idx.replace("rs", "") for idx in ids])
    post = Entrez.read(Entrez.epost(db="snp", id=id_str))
    summaries = Entrez.read(Entrez.esummary(
        db="snp", 
        query_key=post["QueryKey"], 
        webenv=post["WebEnv"], 
        retmax=10000
    ))
    genes = []
    for s in summaries["DocumentSummarySet"]["DocumentSummary"]:
        if "GENES" in s and len(s["GENES"]) > 0:
            genes.append({
                "snp": "rs" + s["SNP_ID"],
                "genes": ",".join(g["NAME"] for g in s["GENES"]),
                "clinical": s["CLINICAL_SIGNIFICANCE"],
                "chromosome": s["CHR"],
                "function_class": s["FXN_CLASS"]
            })
        else:
            genes.append(dict())
    return pd.DataFrame.from_records(genes)


def summarize_associations(df, name):
    """Summarize results for a subset of metaoblites."""
    subset = df.copy()
    n_snvs = subset.variant.nunique()
    n_microbes = subset.Genus.nunique()
    n_metab = subset.metabolite.nunique()

    pl = (
        ggplot(subset, aes(y="r2", x="metabolite", fill="type"))
        + geom_bar(stat="identity") 
        + scale_x_discrete(limits=subset.metabolite[::-1].drop_duplicates())
        + scale_y_continuous(labels=percent_format())
        + coord_flip()
        + labs(y = "explained metabolite variance", x="")
        + scale_fill_manual(values={"genetics_r_squared": "steelblue", "micro_r_squared": "mediumseagreen"})
        + guides(fill = None)
        + theme_minimal() 
        + theme(figure_size=(3, 0.4*n_metab))
    )
    pl.save(f"../figures/{name}_r2.pdf")

    subset_ids = joint_r_sq[joint_r_sq.metabolite.str.contains("subsetchol")].index
    subset_snvs = sig_metab_assoc[sig_metab_assoc.Phenotype.isin(subset_ids)]
    subset_snvs["variant"] = subset_snvs["rsid"] + " | chr" + subset_snvs["chromosome"] 
    subset_snvs.loc[subset_snvs.genes != "", "variant"] = (
        subset_snvs.loc[subset_snvs.genes != "", "variant"] 
        + " | " 
        + subset_snvs.loc[subset_snvs.genes != "", "genes"]
    )
    subset_snvs["metabolite"] = subset_snvs["metabolite"].str.replace(" (1)", "", regex=False)
    subset_snvs["value"] = 1

    subset_mat = subset_snvs.pivot_table(index="metabolite", columns="variant", values="value", fill_value=0)
    snv_map = sns.clustermap(subset_mat, cmap="Blues", figsize=(0.2*n_snvs, 0.2*n_metab), yticklabels=True)
    plt.savefig(f"../figures/{name}_variants.pdf", pad_inches=0.1, bbox_inches="tight")
    
    subset_mic = micro_metab_assoc[micro_metab_assoc.metabolite_id.isin(subset_ids)]
    subset_mic["metabolite"] = subset_mic["metabolite"].str.replace(" (1)", "", regex=False)

    subset_mic_mat = subset_mic.pivot_table(columns="Genus", index="metabolite", values="r", fill_value=0)
    microbe_map = sns.clustermap(subset_mic_mat, cmap="seismic", figsize=(0.2*n_microbes, 0.2*n_metab), yticklabels=True, xticklabels=True, center=0)
    plt.savefig(f"../figures/{name}_genera.pdf", pad_inches=0.1, bbox_inches="tight")
    
    return pl, snv_map, microbe_map