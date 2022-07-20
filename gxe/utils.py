"""Some general utility functions."""

from Bio import Entrez
import matplotlib.pyplot as plt
import pandas as pd
import re
import seaborn as sns
from plotnine import *
from mizani.formatters import percent_format
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import numpy as np

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


def summarize_associations(name, joint_r_sq, sig_metab_assoc, micro_metab_assoc, only_r2=False):
    """Summarize results for a subset of metabolites."""
    subset = (
        joint_r_sq.copy()[["metabolite", "BIOCHEMICAL_NAME", "group", "micro_r2", "geno_r2"]]
        .melt(id_vars=["metabolite", "BIOCHEMICAL_NAME", "group"], value_name="r2", var_name = "type")
        .sort_values(["group", "r2"])
    )
    subset = subset[subset.r2 > 0.0001]
    subset.BIOCHEMICAL_NAME = subset.BIOCHEMICAL_NAME.str.replace(" (1)", "", regex=False)
    n_metab = subset.metabolite.nunique()

    pl = (
        ggplot(subset, aes(y="r2", x="BIOCHEMICAL_NAME", fill="type"))
        + geom_bar(stat="identity") 
        + scale_x_discrete(limits=subset.BIOCHEMICAL_NAME[::-1].drop_duplicates())
        + scale_y_continuous(labels=percent_format())
        + coord_flip()
        + labs(y = "explained metabolite variance", x="")
        + scale_fill_manual(values={"geno_r2": "steelblue", "micro_r2": "mediumseagreen"})
        + guides(fill = None)
        + theme_minimal() 
        + theme(figure_size=(3, 0.25*n_metab))
    )
    pl.save(f"figures/{name}_r2.pdf")
    
    if only_r2:
        return pl

    subset_snvs = sig_metab_assoc[sig_metab_assoc.metabolite.isin(subset.metabolite)]
    subset_snvs["variant"] = subset_snvs["rsid"] + " | chr" + subset_snvs["CHR"].astype(str) 
    subset_snvs.loc[~subset_snvs.genes.isnull(), "variant"] = (
        subset_snvs.loc[~subset_snvs.genes.isnull(), "variant"] 
        + " | " 
        + subset_snvs.loc[~subset_snvs.genes.isnull(), "genes"]
    )
    subset_snvs = pd.merge(subset_snvs, subset[["metabolite", "BIOCHEMICAL_NAME"]].drop_duplicates(), on="metabolite")
    n_snvs = subset_snvs.variant.nunique()

    subset_mat = subset_snvs.pivot_table(index="BIOCHEMICAL_NAME", columns="variant", values="BETA", fill_value=0)
    snv_map = sns.clustermap(
        subset_mat, 
        cmap="seismic", 
        center=0,
        figsize=(6 + 0.15*n_snvs, 3 + 0.15*subset_snvs.metabolite.nunique()), 
        yticklabels=True, 
        xticklabels=True, 
        metric="jaccard",
        row_cluster=subset_mat.shape[0] > 1)
    plt.setp(snv_map.ax_heatmap.get_xticklabels(), rotation=45, ha='right') 
    plt.savefig(f"figures/{name}_variants.pdf", pad_inches=0.1, bbox_inches="tight")
    
    subset_mic = micro_metab_assoc[micro_metab_assoc.metabolite.isin(subset.metabolite)]
    subset_mic = pd.merge(subset_mic, subset[["metabolite", "BIOCHEMICAL_NAME"]].drop_duplicates(), on="metabolite")
    subset_mic["genus"] = subset_mic.taxon.str.split("|").str[1]
    n_microbes = subset_mic.genus.nunique()

    subset_mic_mat = subset_mic.pivot_table(columns="genus", index="BIOCHEMICAL_NAME", values="r", fill_value=0)
    microbe_map = sns.clustermap(subset_mic_mat, cmap="seismic", figsize=(10 + 0.12*n_microbes, 3 + 0.2*subset_mic.BIOCHEMICAL_NAME.nunique()), yticklabels=True, xticklabels=True, center=0)
    plt.setp(microbe_map.ax_heatmap.get_xticklabels(), rotation=45, ha='right') 
    plt.savefig(f"figures/{name}_genera.pdf", pad_inches=0.1, bbox_inches="tight")
    
    return pl, snv_map, microbe_map


def stars(p):
    if p>0.05:
        return "n.s."
    else:
        return sum(p < c for c in [0.05, 0.01, 0.001]) * "٭"

    
def enrichment(data, q_cutoff, column, figsize=(3, 4), min_sig=1):
    """Plot an enrichment analysis for the tests."""
    sig_data = data[data.q < q_cutoff]
    full_counts = data[column].value_counts()
    sig_counts = sig_data[column].value_counts()
    stats = pd.DataFrame({"full_counts": full_counts, "sig_counts": sig_counts, "n": full_counts.sum(), "sig_n": sig_counts.sum()}).fillna(0)
    stats["p"] = stats.apply(lambda df: hypergeom.sf(max(0, df.sig_counts - 1), df.n, max(df.full_counts, 1), df.sig_n), axis=1)
    stats["odds"] = (stats.sig_counts / stats.sig_n) / (stats.full_counts / stats.n)
    stats["log_odds"] = np.log(stats.odds)
    stats["q"] = multipletests(stats.p, method="fdr_bh")[1]
    
    stats["all"] = stats["full_counts"] / stats["n"]
    stats["significant"] = stats["sig_counts"] / stats["sig_n"]
    stats.index.name = column
    stats.reset_index(inplace=True)
    
    long = stats[stats.sig_counts >= min_sig].melt(id_vars=[column, "p", "q"], value_vars=["all", "significant"], value_name="prevalence", var_name="group")
    long[column] = pd.Categorical(long[column], long.sort_values(by="prevalence")[column].unique()) 
    long["sig_stars"] = long.q.apply(stars)
    
    pl = (
        ggplot(long, aes(x="prevalence", y=column, shape="group", color="group")) +
        geom_line(aes(group=column), color="black") +
        geom_point(size=2) +
        theme_minimal() +
        theme(figure_size=figsize) +
        labs(y="")
    )
    if (long.q < 0.05).any():
        pl += geom_text(
            aes(label="sig_stars"), 
            data=long[(long.q<0.05) & (long.group == "significant")], 
            color="black", va="center", ha="left", 
            nudge_x=(long.prevalence.max() - long.prevalence.min())*0.025)
        pl += xlim(0, long.prevalence.max()*1.1)
    
    return stats.sort_values(by="p"), pl