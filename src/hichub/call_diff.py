import pandas as pd
import numpy as np
import igraph as ig
import pybedtools
from scipy import stats
from optparse import OptionParser
import sys, os, multiprocessing
import random
import copy
from statsmodels.stats.multitest import multipletests

random.seed(5)
def Dis_df_hic(_df_interaction,  _col_fore, _col_back, _resolution):
    col_fore = _col_fore
    col_back  = _col_back
    resolution = _resolution
    _df_interaction = _df_interaction.fillna(0)
    

    col1 = _col_fore
    col2 = _col_back
    hic_merge = _df_interaction
    hic_merge['distance'] = (hic_merge['bin2'] - hic_merge['bin1'])/10**4        
    hic_norm = pd.DataFrame(columns=['#chr','bin1','bin2',col1,col2])
    hic_norm_group = hic_merge.groupby('distance')
    for hic_norm_subgroup in hic_norm_group:
        hic_norm_test = hic_norm_subgroup[1]
        hic_norm_test = Distance_Loess_Norm(hic_norm_test, col1, col2)
        hic_norm = pd.concat([hic_norm, hic_norm_test])
            
    df_interaction = hic_norm #LOESS_Norm_df(_df_interaction, col_fore, col_back)
    df_interaction.loc[:,'#chr']=df_interaction.iloc[:,0].astype(str).replace('chr','')
    df_interaction.loc[:,'#chr1']=df_interaction.iloc[:,0]
    df_interaction.loc[:,'x1']=df_interaction.iloc[:,1].astype(int)
    df_interaction.loc[:,'x2']=df_interaction.iloc[:,1].astype(int)+int(resolution)
    df_interaction.loc[:,'chr2']=df_interaction.iloc[:,0]
    df_interaction.loc[:,'y1']=df_interaction.iloc[:,2].astype(int)
    df_interaction.loc[:,'y2']=df_interaction.iloc[:,2].astype(int)+int(resolution)
        
    if('logFC' in df_interaction.columns):
        df_interaction.loc[:,'log_FC'] = df_interaction.loc[:,'logFC']
    else:
        df_interaction.loc[:,'log_FC'] = df_interaction.loc[:,col_fore] - df_interaction.loc[:,col_back]
    df_interaction = df_interaction.loc[:,['#chr1','x1','x2','chr2','y1','y2','log_FC', col_fore, col_back]]
    
    return df_interaction

def Distance_Loess_Norm(df_test, _col1, _col2):
    distance = df_test['distance'].unique()[0]
    
    if(distance<=300):
        df_test[_col1] = df_test[_col1] + 1
        df_test[_col2] = df_test[_col2] + 1
        n_bins = 100
        df_test['A'] = 0.5*(np.log2(df_test[_col1]) + np.log2(df_test[_col2]))
        A_min = df_test['A'].min()
        A_max = df_test['A'].max()
        df_out = pd.DataFrame()
        
        if A_min != A_max :
            x = np.arange(A_min, A_max, (A_max-A_min)/n_bins)
            for i in range(len(x)-1):
                df_bbb = df_test[df_test['A']>=x[i]]
                df_bbb = df_bbb[df_bbb['A']<=x[i+1]]
                sum_1 = df_bbb[_col1].sum(axis=0)
                sum_2 = df_bbb[_col2].sum(axis=0)
                if(sum_1>0 and sum_2>0):
                    df_bbb[_col2] = df_bbb[_col2]/sum_2*sum_1
                df_out = pd.concat([df_out,df_bbb])
            return df_out.loc[:,['#chr','bin1','bin2', _col1, _col2]]
        else:
            sum_1 = df_test[_col1].sum(axis=0)
            sum_2 = df_test[_col2].sum(axis=0)
            df_test[_col2] = df_test[_col2]/sum_2*sum_1
            return df_test.loc[:,['#chr','bin1','bin2', _col1, _col2]]
        
    if(distance>300):
        sum_1 = df_test[_col1].sum(axis=0)
        sum_2 = df_test[_col2].sum(axis=0)
        df_test[_col2] = df_test[_col2]/sum_2*sum_1
        return df_test.loc[:,['#chr','bin1','bin2', _col1, _col2]]
    

def Norm_df_hic(_df_interaction,  _col_fore, _col_back, _resolution):
    col_fore = _col_fore
    col_back  = _col_back
    resolution = _resolution
    _df_interaction = _df_interaction.fillna(0)
    df_interaction = LOESS_Norm_df(_df_interaction, col_fore, col_back)
    df_interaction.loc[:,'#chr']=df_interaction.iloc[:,0].astype(str).replace('chr','')
    df_interaction.loc[:,'#chr1']=df_interaction.iloc[:,0]
    df_interaction.loc[:,'x1']=df_interaction.iloc[:,1].astype(int)
    df_interaction.loc[:,'x2']=df_interaction.iloc[:,1].astype(int)+int(resolution)
    df_interaction.loc[:,'chr2']=df_interaction.iloc[:,0]
    df_interaction.loc[:,'y1']=df_interaction.iloc[:,2].astype(int)
    df_interaction.loc[:,'y2']=df_interaction.iloc[:,2].astype(int)+int(resolution)
        
    if('logFC' in df_interaction.columns):
        df_interaction.loc[:,'log_FC'] = df_interaction.loc[:,'logFC']
    else:
        df_interaction.loc[:,'log_FC'] = df_interaction.loc[:,col_fore] - df_interaction.loc[:,col_back]
    df_interaction = df_interaction.loc[:,['#chr1','x1','x2','chr2','y1','y2','log_FC', col_fore, col_back]]
    
    return df_interaction

def Add_pesudo_count(_df_test, _count, _col1, _col2):
    
    df_test = _df_test
    count = _count
    col1 = _col1
    col2 = _col2
    
    df_test[col1] = df_test[col1] + count
    df_test[col2] = df_test[col2] + count
    
    return df_test

def LOESS_Norm_df (_df, _col1, _col2):
    ## this is a similar approach as LOESS Normalization
    df_test = _df
    df_test = Add_pesudo_count(df_test, 1, _col1, _col2)
    n_bins = 100
    df_test['A'] = 0.5*(np.log2(df_test[_col1]) + np.log2(df_test[_col2])) ## A is the value for MA plot 
    
    A_min = df_test['A'].min()
    A_max = df_test['A'].max()

    x = np.arange(A_min, A_max, (A_max-A_min)/n_bins)

    df_out = pd.DataFrame()
    for i in range(len(x)-1):
        df_bbb = df_test[df_test['A']>=x[i]]
        df_bbb = df_bbb[df_bbb['A']<x[i+1]]
        sum_1 = df_bbb[_col1].sum(axis=0)
        sum_2 = df_bbb[_col2].sum(axis=0)
        df_bbb[_col2] = round(df_bbb[_col2]/sum_2*sum_1, 2)
       
        df_out = pd.concat([df_out,df_bbb],axis=0)#df_out.append(df_bbb)
        df_out = df_out.sort_index()
        
    return df_out.sort_index()

def Convert_Loops_to_Graph(_df_hic, _weight_col):
	## Assign a list of weight ot graph
	
	## loop format: ['#chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'GeneID', 'weight_cols']
	df_bins = Loops_Return_two_bins_no_dup(_df_hic)
	
	## eliminate float in chr
	df_bins['name'] = df_bins['#chr1'].astype(str).str.split(".",expand=True)[0]+':'+df_bins['x1'].astype(int).astype(str)+'-'+df_bins['x2'].astype(int).astype(str)
	Num_vs = len(df_bins.index)
	## Initiation a graph from loops file 
	graph_tem = ig.Graph()
	graph_tem.add_vertices(Num_vs)
	graph_tem.vs["name"] = df_bins.loc[:,'name']
	df_edge = _df_hic.merge(df_bins, on=['#chr1', 'x1', 'x2']).merge(
		df_bins, left_on=['chr2', 'y1', 'y2'], right_on=['#chr1', 'x1', 'x2'])
	graph_tem.add_edges(df_edge.loc[:, ['index_x','index_y']].values)

	for weight in _weight_col:
		if (weight in _df_hic.columns):
			graph_tem.es[weight] = df_edge.loc[:,weight].values
            
	return graph_tem

def Loops_Return_two_bins_no_dup(df_hic):
	## Associated by promoter
	second_bin_columns = [3,4,5,0,1,2]+list(range(6,len(df_hic.columns),1))
	df_hic=pd.concat([df_hic,pd.DataFrame(df_hic.iloc[:, second_bin_columns].values, columns=df_hic.columns)],axis=0).sort_index()#df_hic.append(pd.DataFrame(df_hic.iloc[:, second_bin_columns].values, columns=df_hic.columns),sort=False).sort_index()
	return df_hic.iloc[:,0:3].drop_duplicates().reset_index().drop('index',axis=1).reset_index()

def convert_cluster2bed(df_cluster, usecol):
	df_tem = df_cluster[usecol].str.split(r"\:|-",expand=True)
	df_tem = pd.concat( [df_tem, df_cluster], axis=1)
	if (df_tem.iloc[0,0].find('chr') == -1):
		df_tem[0] = 'chr'+df_tem[0]
        
	return df_tem

def convert_bin2bed(df_cluster, col_name):
	df_tem = df_cluster[col_name].str.split(r"\:|-",expand=True)
	df_tem = pd.concat( [df_tem, df_cluster], axis=1)
	if (df_tem.iloc[0,0].find('chr') == -1):
		df_tem[0] = 'chr'+df_tem[0]
        
	return df_tem

def convert_vs2bed(input_graph, col_name):
	## output first 3 columns is standard bed format
	df_tem = pd.DataFrame(data={col_name:input_graph.vs[col_name]})
	df_tem = pd.concat( [df_tem[col_name].str.split(r"\:|-",expand=True),df_tem], axis=1)
	if (df_tem.iloc[0,0].find('chr') == -1):
		df_tem[0] = 'chr'+df_tem[0]
        
	return df_tem

def convert_graph_vs_to_df(_input_graph):
	df_vs = pd.DataFrame(data= {"degree":_input_graph.degree()})
	for col in _input_graph.vs.attributes():
		df_vs[col] = _input_graph.vs[col]
        
	return df_vs

def Graph_Pagerank(_input_graph):
	input_graph = _input_graph
	input_graph.vs['pagerank'] = input_graph.pagerank(weights=input_graph.es['weight'])
    
	return input_graph

def Stich_Region_Above_global_Mean(_graph, _resolution, _gap_size, _mean):
	resolution=_resolution
	df_vs_graph = convert_graph_vs_to_df(_graph)
	df_nodes = convert_cluster2bed(df_vs_graph, 'name')
	df_nodes[1] = df_nodes[1].astype(int)
	df_nodes = df_nodes.sort_values(by=1)
	df_nodes = df_nodes[df_nodes['pagerank'] > _mean] ## Only use nodes > mean
	Report_list=[]
	if (len(df_nodes)>0):
		## report stich regions
		
		reg_chr = str(df_nodes.iloc[0,0])
		reg_start= int(df_nodes.iloc[0,1])
		reg_end = int(reg_start)

		for bin1 in df_nodes.iloc[:,1].astype(int):
			if (bin1-reg_end)<=_gap_size*resolution:
				reg_end = bin1
			else:
				Report_list.append([reg_chr+':'+str(reg_start)+'-'+str(reg_end), _gap_size])
				reg_start = bin1
				reg_end = bin1
		Report_list.append([reg_chr+':'+str(reg_start)+'-'+str(reg_end), _gap_size])
        
	return pd.DataFrame(data=Report_list, columns=['hub_name', 'merge_level'])

def Return_Sorted_Adjacency_Matrix(_graph, _attr):
	graph_tem = _graph
	attr      =_attr
	idx_name = [int(str(x).split(":")[1].split("-")[0]) for x in graph_tem.vs['name']]
	
	matrix_tem = pd.DataFrame(data=graph_tem.get_adjacency(attribute=attr), columns=idx_name, index=idx_name)
	df_reindex = pd.DataFrame(data={ 'rank': (stats.rankdata(matrix_tem.columns)-1).astype(int)})
	idx_rank = df_reindex.sort_values(by='rank').index
    
	return matrix_tem.iloc[idx_rank, :].T.iloc[idx_rank, :]

def Pvalue_Rank_Test_Matrix(_matirx):
    matrix_for_test = _matirx
    data_test = matrix_for_test.values.flatten()
    data_test = data_test[~ np.isnan(data_test)]
    if (len(data_test)>10):
        w, pvalue =stats.wilcoxon(data_test, zero_method='zsplit', alternative='greater', correction=True, mode='approx')  
        # “zsplit”: Includes zero-differences in the ranking process and split the zero rank between positive and negative ones.
    else:
        pvalue=1.0
        
    return float(pvalue)

def Return_Pvalue_For_Given_Graph(_df_region, _resolution, _matrix):
    df_region = _df_region
    df_regionh_bed = convert_cluster2bed(df_region, 'hub_name').sort_values(by=1)
    resolution = _resolution
    matrix_for_test = _matrix
    
    idx_regs = []
    for name_stitch in df_region.hub_name:
        region_loc= name_stitch.split(":")[1].split("-")
        idx_reg = []
        for idx in matrix_for_test.index:
            if ((idx>=int(region_loc[0]))&(idx<=int(region_loc[1]))):
                idx_reg.append(idx)
        idx_regs.append(idx_reg)

    pvalue_region= []
    for i in range(len(idx_regs)):
        part_matrix_for_test = matrix_for_test.loc[idx_regs[i],:].T.loc[idx_regs[i], :]
        pvalue_tem = Pvalue_Rank_Test_Matrix(part_matrix_for_test)
        
        if (pvalue_tem < 1):
            if (pvalue_tem == 0):
                pvalue_tem = 1e-400
            pvalue_region.append([df_region.hub_name[i],df_region.hub_name[i],-np.log10(pvalue_tem)])
            for j in range(0, i, 1):
                part_matrix_for_test = matrix_for_test.loc[idx_regs[i],:].T.loc[idx_regs[j], :]
                pvalue_tem = Pvalue_Rank_Test_Matrix(part_matrix_for_test)
                pvalue_region.append([df_region.hub_name[i],df_region.hub_name[j], np.round(-np.log10(pvalue_tem),3)])

    return pd.DataFrame(data=pvalue_region, columns=['left_hub_anchor', 'right_hub_anchor', '-log10(pvalue)']).sort_values('-log10(pvalue)', ascending=False)

def Main_For_Diff_Regions(df_hic, _col_fore, _col_back,  _resolution, _pvalue, _logFC):
    data_1 = copy.deepcopy(df_hic)
    data_ori = revise_data_format(data_1, _col_fore, _col_back)    
    _gapsize=2
    logfc_cutoff=_logFC
    cut_pvalue=-np.log10(_pvalue)
    _df_hic = df_hic
    _df_hic[_col_fore+'_weight'] = _df_hic[_col_fore]*_df_hic.log_FC.apply(lambda x: 1 if x > logfc_cutoff else(0))  
    Norm_window_Size=0    
    _df_hic['diff'] = _df_hic[_col_fore] - _df_hic[_col_back]
    _df_hic['pagerank_weight'] = _df_hic['diff']*(abs(_df_hic.y1-_df_hic.x1)).apply(lambda x : 1 if x >= Norm_window_Size*_resolution else (0) )
    _df_hic['pagerank_weight'] = _df_hic['pagerank_weight'].apply(lambda x : x if x >0 else (0) )    
    weight_list= ['diff','pagerank_weight', _col_fore+'_weight']
    input_graph = Convert_Loops_to_Graph(_df_hic, weight_list)   
    input_graph.es['weight'] = input_graph.es['pagerank_weight']
    input_graph = Graph_Pagerank(input_graph)
    global_median = np.percentile(input_graph.vs['pagerank'], 50)
    cut_off = global_median
    input_graph.es['weight'] = input_graph.es[_col_fore+'_weight']
    structure = input_graph.community_multilevel(weights=input_graph.es['weight'], return_levels=True)
    df_out = pd.DataFrame(columns=['left_hub_anchor', 'right_hub_anchor', '-log10(pvalue)'])
    
    if (len(structure)==0):
        df_out.to_csv(_col_back+'_'+_col_fore+'_specific_regions.bed', sep='\t', mode='a', header=False, index=None)
        return None
    i=0
    sub_structure = structure[0].subgraphs()
    for graph_tem in sub_structure:
        if (len(graph_tem.vs)>_gapsize+1):
            i+=1
            df_hubs = Stich_Region_Above_global_Mean(graph_tem, _resolution, _gapsize, cut_off)
            if(len(df_hubs)>0):
                Diff_matrix = Return_Sorted_Adjacency_Matrix(graph_tem, 'diff')
                Diff_matrix = Diff_matrix.where(np.tril(np.ones(Diff_matrix.shape)).astype(bool))
                df_out = pd.concat([df_out,Return_Pvalue_For_Given_Graph(df_hubs, _resolution, Diff_matrix)],axis=0)#df_out.append(Return_Pvalue_For_Given_Graph(df_hubs, _resolution, Diff_matrix))
                
    #data_ori = revise_data_format(df_hic, _col_fore, _col_back)
    df_test = copy.deepcopy(df_out)
    p_value_new = re_calculate_pvalue(df_test, data_ori)
    #print(p_value_new)
    df_out.loc[:,'-log10(pvalue)_new'] = np.round(-np.log10(p_value_new),3)
    #print(df_out)
    df_out['-log10(pvalue)'] = np.round(df_out['-log10(pvalue)'].astype(float),3)
    df_out = df_out.sort_values(by='-log10(pvalue)', ascending=False)
    df_out = df_out[df_out['-log10(pvalue)_new']>cut_pvalue]
    df_out.to_csv(_col_back+'_'+_col_fore+'_specific_regions.bed', sep='\t', mode='a', header=False, index=None)
    
    return structure

def revise_data_format(df_hic, _col_fore, _col_back):
    data_ori = df_hic
    data_ori['diff'] = data_ori[_col_fore] - data_ori[_col_back]
    data_ori.loc[:,'#chr'] = data_ori.iloc[:,0].astype(str).replace('chr','')
    data_ori = data_ori.loc[:,['#chr1','x1','y1','log_FC','diff']]
    
    return data_ori

def re_calculate_pvalue(_df,data):
    if len(_df) == 0:
        return []
    test = _df
    test['chr_x'] = test['left_hub_anchor'].str.split(':', expand = True)[0]
    test['bin1_x'] = test['left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    test['bin2_x'] = test['left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    test['chr_y'] = test['right_hub_anchor'].str.split(':', expand = True)[0]
    test['bin1_y'] = test['right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    test['bin2_y'] = test['right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    test = test.loc[:,['chr_x','bin1_x','bin2_x','chr_y','bin1_y','bin2_y']]
    p_value_new = count_number(test, data)
    
    return p_value_new
    
def count_number(df1,test):
    pvalue_box = []
    region = df1
    for i in range(len(df1)):
        data = test
        if region.iloc[i,1]<=region.iloc[i,4]:        
            data = data[data['x1']>=region.iloc[i,1]]
            data = data[data['x1']<=(region.iloc[i,2])]
            data = data[data['y1']>=region.iloc[i,4]]
            data = data[data['y1']<=(region.iloc[i,5])]
        else:
            data = data[data['x1']>=region.iloc[i,4]]
            data = data[data['x1']<=(region.iloc[i,5])]
            data = data[data['y1']>=region.iloc[i,1]]
            data = data[data['y1']<=(region.iloc[i,2])]
            
        if len(data)>10:           
            w, pvalue =stats.wilcoxon(data['diff'], zero_method='zsplit', alternative='greater', correction=True, mode='approx')
            pvalue_box.append(pvalue)
        else:
            pvalue_box.append(1)
        
        
    return pvalue_box

def multi_task(_chr_name, _df_chr, _col_fore, _col_back,  _resolution, _pvalue, _logFC, norm):
    col_fore = _col_fore
    col_back = _col_back
    resolution = _resolution
    pvalue = _pvalue
    logFC = _logFC
    if norm == 'LOESS':
        df_chr = Norm_df_hic(_df_chr, col_fore, col_back, resolution)
    elif norm == 'Dis_LOESS':
        df_chr = Dis_df_hic(_df_chr, col_fore, col_back, resolution)
    else:
        print('Please enter one of the following two normalization method:')
        print('1) LOESS')
        print('2) Dis_LOESS')
    structure = Main_For_Diff_Regions(df_chr, col_fore, col_back, resolution, pvalue, logFC)
    
    return structure

def cut_off_value(_df_test, _value, _col1, _col2):
    
    df_test = _df_test
    value = _value
    col1 = _col1
    col2 = _col2
    
    col_list = df_test.columns.tolist()
    
    df_test[col1 + r'+' + col2] = df_test[col1] + df_test[col2]
    df_test = df_test[(df_test[col1+r'+'+col2] > value)]
    
    df_test = df_test.loc[:,[col_list[0], col_list[1], col_list[2], col1, col2]]
    
    return df_test

def Multi_Main_For_Diff_Regions(_PATH_interaction, _col_fore, _col_back,  _resolution, _pvalue, _num_threads, _logFC, _cut_off, norm):
    if True:
        PATH_interaction = _PATH_interaction
        
        test = pd.read_csv(PATH_interaction, sep='\t', dtype={'#chr':str}, iterator=True, chunksize=100000)


        out = pd.DataFrame()
        for chunk in test:
            chunk = chunk.fillna(0)
            chunk = chunk[chunk[_col_fore]+chunk[_col_back] > _cut_off]
            out = pd.concat([out,chunk],axis=0)
        
        df_tem = out
        
        #df_tem = df_tem[df_tem['#chr'] == '1']

        logFC = _logFC
        col_fore=_col_fore
        col_back=_col_back
        resolution = _resolution
        _col = df_tem.columns[0]
        df_groups = df_tem.groupby(by=_col)
        
        for df_group in df_groups:
            chr_name = df_group[0]
            df_hic_chr = df_group[1]
            structure = multi_task(chr_name, df_hic_chr, col_fore, col_back, resolution, _pvalue, logFC, norm)
            
            i = 0
            if structure is not None:
                for graph_tem in structure[0].subgraphs():
                    df_test = convert_graph_vs_to_df(graph_tem)
                    df_test['cluster'] = 'cluster'+str(i)
                    df_test['node_chr'] = df_test['name'].str.split(':', expand = True)[0]
                    df_test['node_p'] = df_test['name'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
                    df_test = df_test.loc[:,['node_chr','node_p','cluster']]
                    df_test.to_csv('cluster_'+ col_fore +'.txt',sep='\t',index=None,mode='a', header = False)
                    i=i+1
            else:
                df_test = pd.DataFrame(data = {'node_chr':[],'node_p':[],'cluster':[]}).to_csv('cluster_'+ col_fore +'.txt',sep='\t',index=None,mode='a', header = False)
            
        print('All subprocesses done.')

        col_out = [ 'left_hub_anchor', 'right_hub_anchor' ,'-log10(pvalue)','-log10(pvalue)_new']
        df_output = pd.read_csv(col_back+'_'+col_fore+'_specific_regions.bed', sep='\t', header=None, names=col_out)
        df_output.sort_values(by='-log10(pvalue)_new', ascending=False).to_csv(str(len(df_output))+'_'+str(logFC)+'_'+col_back+'_'+col_fore+'_specific_regions.bed', sep='\t', mode='a', header=True, index=None)
        
        df_out = pd.read_csv(str(len(df_output))+'_'+str(logFC)+'_'+col_back+'_'+col_fore+'_specific_regions.bed', sep='\t')
        df_out_p = pow(10, -df_out['-log10(pvalue)_new'])

        if len(df_out) == 0:
            df_out_p_new = []
            df_out['-log10(FDR)'] = []
        else:
            df_out_p_new = multipletests(df_out_p, method='fdr_bh', is_sorted=False, returnsorted=False)
            df_out['-log10(FDR)'] = -np.log10(df_out_p_new[1])
            
        df_out['-log10(pvalue)'] = df_out['-log10(pvalue)_new']
        df_out = df_out.loc[:,['left_hub_anchor','right_hub_anchor','-log10(pvalue)','-log10(FDR)']]
        df_out['haha'] = df_out['left_hub_anchor']
        df_out['left_hub_anchor'] = df_out['right_hub_anchor']
        df_out['right_hub_anchor'] = df_out['haha']
        df_out = df_out.loc[:,['left_hub_anchor','right_hub_anchor','-log10(pvalue)','-log10(FDR)']]
        os.remove(str(len(df_output))+'_'+str(logFC)+'_'+col_back+'_'+col_fore+'_specific_regions.bed')
        df_out.to_csv(col_fore+'_specific_hubs.bed', sep='\t', index=None, header=True)
        
        os.remove(col_back+'_'+col_fore+'_specific_regions.bed')
        
    return None

def revise_output_format(state,name):
    input_file = pd.read_csv(state+'_specific_hubs.bed', sep='\t')
    
    input_file['#chr'] = input_file['left_hub_anchor'].str.split(':', expand = True)[0]
    input_file['left_hub_anchor_start'] = input_file['left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    input_file['left_hub_anchor_end'] = input_file['left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)
    input_file['right_hub_anchor_start'] = input_file['right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    input_file['right_hub_anchor_end'] = input_file['right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)
    input_file = input_file.loc[:,['#chr', 'left_hub_anchor_start', 'left_hub_anchor_end', 'right_hub_anchor_start', 'right_hub_anchor_end', '-log10(pvalue)','-log10(FDR)']]
    input_file.to_csv(state+'_specific_hubs_comparing_with_'+name+'.bed', sep='\t', index=None, header=True)
    os.remove(state+'_specific_hubs.bed')
    return None

### End of Visulization
####################################################################################
### FUNCTION
### FUNCTIONS
def run(argv):
	## parameters
	PATH_INPUT=argv.input_path
	
	col_fore = argv.label.split(',')[0]
	col_back  = argv.label.split(',')[1]
	
	resolution = argv.res
	pvalue=argv.pvalue
	logFC = argv.foldchange
	cut_off = argv.cut_off_value
	num_threads=argv.thread
	norm = argv.norm_method
	
	print (" ")
	#print("Run main")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in txt format: %s" % PATH_INPUT)
	print ("First label is: %s" % col_fore)
	print ("Second label is: %s" % col_back)
	print ("Resolution is: %i" % resolution)
	print ("Normalizaiton method is: %s" % norm)
	print ("Pvalue cutoff for output (diff hub) is: %s" % pvalue)
	print ("FC for qualified difference is: %s" % logFC)
	print ("Threshold value for sum of row counts is: %s" % cut_off)
	#print ("Number of threads used is: %i" % num_threads)
	print ("End of Summary.")
	print (" ")
	
#### Main 
	Multi_Main_For_Diff_Regions(PATH_INPUT, col_fore, col_back, resolution, pvalue, num_threads, logFC, cut_off, norm)		
	revise_output_format(col_fore)
	Multi_Main_For_Diff_Regions(PATH_INPUT, col_back, col_fore, resolution, pvalue, num_threads, logFC, cut_off, norm)
	revise_output_format(col_back)


	print(" ")
    
	return None

def main(argv):
	desc="Main program to call hubs."
	
	parser = OptionParser(description=desc)
	
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Input .txt file you want to call hubs.", metavar="<file>")
	
	parser.add_option("-l", "--label", action="store", type="string",
			dest="label", help="Name of labels you want to call hubs in your .txt files. Delimiter ',', for example 'label1,label2'.", metavar="<str>")

	parser.add_option("-r", "--resolution", action="store", type="int",
		dest="res", help="Resolution of .hic files you converted before.", metavar="<int>")  
      
	parser.add_option("-d", "--diff", action="store", type="float", default = 0,
		dest="foldchange", help="Optional: difference to evaluate qualified difference, default=0.", metavar="<float>")
	parser.add_option("-c", "--cut_off", action="store", type="float", default = 10,
		dest="cut_off_value", help="Optional: Remove sum of row counts smaller than this value, default=10.", metavar="<float>")
	parser.add_option("-p", "--pvalue", action="store", type="float", default =0.00001,
		dest="pvalue", help="Optional: Pvalue threshold for output (diff hub), default=1e-5.", metavar="<float>")
	
	
	parser.add_option("-t", "--num_threads", action="store", type="int", default =1,
		dest="thread", help="Optional: Number of threads to run, default=1.", metavar="<int>")
    
	parser.add_option("-n", "--normalization_method", action="store", type="str", default = 'LOESS',
		dest="norm_method", help="Normalization method.", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 3:
		parser.print_help()
		sys.exit(1)
	
	## parameters
	PATH_INPUT=opt.input_path
	
	col_fore = opt.label.split(',')[0]
	col_back  = opt.label.split(',')[1]
	
	resolution = opt.res
	pvalue=opt.pvalue
	logFC = opt.foldchange
	cut_off = opt.cut_off_value
	num_threads=opt.thread
	norm = opt.norm_method
    
	print (" ")
	#print("Run main")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in txt format: %s" % PATH_INPUT)
	print ("First label is: %s" % col_fore)
	print ("Second label is: %s" % col_back)
	print ("Resolution is: %i" % resolution)
    
	print ("Normalizaiton method is: %s" % norm)
    
	print ("Pvalue cutoff for output (diff hub) is: %s" % pvalue)
	print ("FC for qualified difference is: %s" % logFC)
	print ("Threshold value for sum of row counts is: %s" % cut_off)
	#print ("Number of threads used is: %i" % num_threads)
	print ("End of Summary.")
	print (" ")
	


#### Main 
	Multi_Main_For_Diff_Regions(PATH_INPUT, col_fore, col_back, resolution, pvalue, num_threads, logFC, cut_off, norm)		
	revise_output_format(col_fore, col_back)
	Multi_Main_For_Diff_Regions(PATH_INPUT, col_back, col_fore, resolution, pvalue, num_threads, logFC, cut_off, norm)
	revise_output_format(col_back, col_fore)

	
	print(" ")
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
