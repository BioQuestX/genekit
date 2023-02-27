import re
import numpy as np
import pandas as pd
from bioinfokit import analys
from bioquest.tl._IO import read_csv_gz
import bioquest as bq

def geneIDconverter(frame, from_id='Ensembl', to_id='Symbol', keep_from=False, gene_type=None):
    annot = read_csv_gz('HumanSymbolEnsemblGencodeV42.csv.gz',index_col=from_id)
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

def unique_exprs(frame,reductions=np.median):
	"""
	基因去重复
	frame: row_index = genenames, columns = samples
	"""
	frame['Ref'] = frame.apply(reductions,axis=1)
	frame.sort_values(by='Ref',ascending=False,inplace=True)
	frame.drop(columns='Ref',inplace=True)
	frame['Ref'] = frame.index
	frame.drop_duplicates(subset='Ref',inplace=True)
	return frame.drop(columns='Ref')
def count2tpm(frame,geneid='Symbol'):
    annot = read_csv_gz('HumanExonLengthGencodeV42.csv.gz',usecols=[geneid,'Length'],index_col=geneid)
    _df = pd.merge(annot, frame, left_index=True, right_index=True)
    nm = analys.norm()
    nm.tpm(df=_df, gl='Length')
    return nm.tpm_norm

def get_TCGA_mRNA(arrow,dtype='tpm',label="01A",gene_id='gene_name',gene_type=None,barcode_length=16):
	_df = pd.read_feather(arrow)
	_df.set_index(keys=gene_id,drop=True,inplace=True)
	# 筛选编码蛋白基因
	if gene_type:
		_df = _df.loc[_df.gene_type==gene_type,:]
	# 筛选数据
	_df = bq.tl.select(frame=_df,pattern=re.compile(pattern=f"^{dtype}_"))

	_df.columns = bq.st.removes(string=_df.columns.values,pattern=re.compile(pattern=f"^{dtype}_"))
    
	_df = bq.tl.select(frame=_df,pattern=re.compile(pattern=f"^.{{13}}{label}"))
    # barcode_length of samples
	_df.columns = bq.st.subs(string=_df.columns.values,start=0,stop=barcode_length)
	# 基因去重复
	return unique_exprs(_df)

# 基因长度https://www.jianshu.com/p/abea4033b61e

