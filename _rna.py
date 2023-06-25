import re
import numpy as np
import pandas as pd
from bioinfokit import analys
import bioquest


def geneIDconverter(frame, from_id='Ensembl', to_id='Symbol',species="Human", keep_from=False, gene_type=None):
    if species == "Human":
        file_path = "h38_gene_info_v43.csv.gz"
    if species == "Mouse":
        file_path = "m39_gene_info_v32.csv.gz"
    annot = bioquest.tl._IO.read_csv_gz(file_path)
    if from_id=='Ensembl':
        annot.loc[:,"Ensembl"] = bioquest.st.removes(string=annot.Ensembl,pattern=r"\.\d+")
    annot.set_index(keys=from_id,inplace=True,drop=True)
    if gene_type:
        gene_type = annot.GeneType.isin(gene_type)
        annot = annot.loc[gene_type:, [to_id]]
    else:
        annot = annot.loc[:, [to_id]]
    if keep_from:
        annoted = pd.merge(annot, frame, left_index=True, right_index=True)
    else:
        annoted = pd.merge(annot, frame, left_index=True, right_index=True)
        annoted.set_index(keys=to_id, inplace=True)
    return annoted



def countto(frame, towhat="tpm", geneid='Ensembl', species="Human"):
    '''
    towhat: tpm(default), fpkm, cpm
    return: a dataframe
    '''
    if species == "Human":
        file_path = "h38_gene_info_v43.csv.gz"
    if species == "Mouse":
        file_path = "m39_gene_info_v32.csv.gz"
    annot = bioquest.tl._IO.read_csv_gz(file_path, usecols=[
                        geneid, 'Length'])
    if geneid=='Ensembl':
        annot.loc[:,"Ensembl"] = bioquest.st.removes(string=annot.Ensembl,pattern=r"\.\d+")
    annot.set_index(keys=geneid,inplace=True,drop=True)
    _df = pd.merge(annot, frame, left_index=True, right_index=True)
    nm = analys.norm()
    nm.tpm(df=_df, gl='Length')
    nm.rpkm(df=_df, gl='Length')
    nm.cpm(df=_df)
    if towhat=="tpm":
        return nm.tpm_norm
    if towhat=="fpkm":
        return nm.rpkm_norm
    if towhat=="cpm":
        return nm.cpm_norm.drop(columns="Length")

