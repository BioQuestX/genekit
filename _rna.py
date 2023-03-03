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


def unique_exprs(frame, reductions=np.median):
    """
    基因去重复
    frame: row_index = genenames, columns = samples
    """
    frame['Ref'] = frame.apply(reductions, axis=1)
    frame.sort_values(by='Ref', ascending=False, inplace=True)
    frame.drop(columns='Ref', inplace=True)
    frame['Ref'] = frame.index
    frame.drop_duplicates(subset='Ref', inplace=True)
    return frame.drop(columns='Ref')


def count2tpm(frame, geneid='Ensembl', species="Human"):
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
    return nm.tpm_norm


def get_TCGA_mRNA(arrow, dtype='tpm', label="01A", gene_id='gene_name', gene_type=None, barcode_length=16):
    _df = pd.read_feather(arrow)
    _df.set_index(keys=gene_id, drop=True, inplace=True)
    # 筛选编码蛋白基因
    if gene_type:
        _df = _df.loc[_df.gene_type == gene_type, :]
    # 筛选数据
    _df = bioquest.tl.select(frame=_df, pattern=re.compile(pattern=f"^{dtype}_"))

    _df.columns = bioquest.st.removes(
        string=_df.columns.values, pattern=re.compile(pattern=f"^{dtype}_"))

    _df = bioquest.tl.select(frame=_df, pattern=re.compile(
        pattern=f"^.{{13}}{label}"))
    # barcode_length of samples
    _df.columns = bioquest.st.subs(string=_df.columns.values,
                             start=0, stop=barcode_length)
    # 基因去重复
    return unique_exprs(_df)

# 基因长度https://www.jianshu.com/p/abea4033b61e
