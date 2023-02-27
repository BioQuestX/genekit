import functools
import pandas as pd
import numpy as np
import plotly.express as px

def _volcano_color(x,data,log2FC=0.585, p = 1, p_adj = 0.05,p_key='pvals_adj',log2FC_key='logfoldchanges'):
	if data[p_key].values[x] < np.min([p,p_adj]):
		if data[log2FC_key].values[x] < -log2FC:
			return 'down'
		elif data[log2FC_key][x] > log2FC:
			return 'up'
		else:
			return 'ns'
	else:
		return 'ns'

def plot_volcano(data,
	p_key='pvals_adj',
	log2FC_key='logfoldchanges',
	p = 1, 
	p_adj = 0.05,
	log2FC=0.585,
	xlabel = 'log2FC',
	ylabel = '-log10(p)',
	od='./'
	):
	_df = data.copy()
	_get_volcano_color = functools.partial(_volcano_color,data=_df,log2FC=log2FC,p=p,p_adj=p_adj,p_key='pvals_adj',log2FC_key='logfoldchanges')
	_df['Color'] = list(map(_get_volcano_color,np.arange(_df.shape[0])))
	_df[p_key] = -np.log10(_df[p_key])
	_df.rename(columns={p_key:ylabel,log2FC_key:xlabel},inplace=True)

	fig = px.scatter(_df, x=xlabel, y=ylabel,color='Color',
		color_discrete_map={'ns':'grey','up':'red','down':'blue'}
	# ,range_y = [0,1000]
	# ,range_x = [-20,20]
	)
	fig.add_hline(y=-np.log10(p_adj),line_dash="dash", line_color="grey")
	fig.add_vline(x=0.585,line_dash="dash", line_color="grey")
	fig.add_vline(x=-0.585,line_dash="dash", line_color="grey")
	fig.write_image(od+'/volcano.pdf',format='pdf',width=800,height=800)
	return fig.show()
