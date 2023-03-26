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
    return stat_res.summary().rename(columns={"log2FoldChange": "logFC"})
