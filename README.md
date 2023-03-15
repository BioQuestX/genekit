Python functions for Bioinformatics

## intall dependency

```bash
conda install -y -c conda-forge numpy matplotlib seaborn pandas matplotlib-venn adjustText tabulate textwrap3
pip install pydeseq2
conda install -y -c bioconda bioinfokit
```

## usage

```python
import genekit as gk
# deferential expression analysis
gk.deseq(count_df,clinical_df,reference,n_jobs)
gk.deg_siglabel(df, lfc='log2FC', lfc_thr=(.585, .585),pv='pvalue',pv_thr=(0.05, 0.05),siglabel=('significant down', 'not significant', 'significant up')) 
gk.deg_filterdeg_filter(df,lfc='log2FC',top_n=None,filter_label='significant down',siglabel=('significant down', 'not significant', 'significant up'))

# convert gene ID to Ensembl/Symbol
gk.geneIDconverter(frame, from_id='Ensembl', to_id='Symbol',species="Human", keep_from=False, gene_type=None)

# remove duplicated genes
gk.unique_exprs(frame, reductions=np.median)

# convert Count to TPM/FPKM/CPM
gk.countto(frame, towhat="tpm", geneid='Ensembl', species="Human")
```