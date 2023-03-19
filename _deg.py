import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import bioquest as bq

def deseq(
        count_df: pd.DataFrame,
        clinical_df: pd.DataFrame,
        reference: str,
        n_jobs: int = 16) -> pd.DataFrame:
    """
    differential expression analysis (DEA) with bulk RNA-seq data
    """
    count_df.index.name = None
    clinical_df.index.name = None
    # 构建DeseqDataSet 对象
    dds = DeseqDataSet(
        counts=count_df,
        clinical=clinical_df,
        reference_level=reference,
        design_factor=clinical_df.columns.values[0],
        refit_cooks=True,
        n_cpus=n_jobs,
    )
    # 离散度和log fold-change评估.
    dds.deseq2()
    # 差异表达统计检验分析
    stat_res = DeseqStats(dds, alpha=0.05, cooks_filter=True,
                          independent_filter=True, n_cpus=n_jobs)
    return stat_res.summary().rename(columns={"log2FoldChange": "log2FC"})


def deg_siglabel(
        df: pd.DataFrame, 
        lfc='log2FC',
        padj:str = None,
        pvalue:str = None,
        lfc_thr=(.585, .585),
        pv_thr=(0.05, 0.05),
        siglabel=('significant down', 'not significant', 'significant up')
        ) -> pd.DataFrame:
    """
    label genes for significant up/down or not significant
    lfc
    pv
    padj
    """
    pv = pvalue if pvalue else padj
    # upregulated
    lg_up = np.logical_and(df[lfc] >= lfc_thr[1],df[pv] < pv_thr[1])
    df.loc[lg_up,'Change'] = siglabel[2] 
    # downregulated
    lg_down = np.logical_and(df[lfc] <= -lfc_thr[0],df[pv] < pv_thr[0])
    df.loc[lg_down, 'Change'] = siglabel[0]
    df.fillna(value={'Change': siglabel[1]}, inplace=True)
    # return df
    print(f'All degs: {df.loc[df.Change != siglabel[1], :].shape[0]}')
    print(
        f'significant up: {df.loc[df.Change == siglabel[2],:].shape[0]}')
    print(
        f'significant down: {df.loc[df.Change == siglabel[0],:].shape[0]}')
    return df


def deg_filter(
        df: pd.DataFrame,
        lfc :str='log2FC',
        top_n=None,
        siglabel=('significant down', 'not significant', 'significant up'),
        filter_label=['significant down']
        ) -> pd.DataFrame:
    if top_n:
        _df = df.sort_values(by=[lfc])
        nrow = df.shape[0]
        dfslice = list(range(0, top_n)) + \
            list(range(nrow-top_n, nrow))
        return _df.iloc[dfslice, :]
    else:
        return bq.tl.subset(df,{"Change":filter_label})
