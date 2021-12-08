"""Some general utility functions."""

from Bio import Entrez
import pandas as pd
import re

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