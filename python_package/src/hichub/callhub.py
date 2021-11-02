#!/usr/bin/env python
########################################################################
## 01/11/2020
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Version.beta
########################################################################
# Usage 
#python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} \
#	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}
########################################################################

import pandas as pd
import numpy as np
import igraph as ig
from pybedtools import BedTool
from scipy import stats
from optparse import OptionParser
import sys
####################################################################################
## FUNCTIONS
### FUNCTION
def convert_loops_to_graph(df_loops, weight_col, _extra_edge_col):
    ## loop format: ['#chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'GeneID', 'weight_cols']
    df_bins = Loops_Return_two_bins_no_dup(df_loops)
    df_bins['name'] = df_bins['#chr1'].astype(int).astype(str)+':'+df_bins['x1'].astype(int).astype(str)+'-'+df_bins['x2'].astype(int).astype(str)
    Num_vs = len(df_bins.index)
    ## Initiation a graph from loops file 
    graph = ig.Graph()
    graph.add_vertices(Num_vs)
    graph.vs["name"] = df_bins.loc[:,'name']
    
    df_edge = df_loops.merge(df_bins, on=['#chr1', 'x1', 'x2']).merge(
        df_bins, left_on=['chr2', 'y1', 'y2'], right_on=['#chr1', 'x1', 'x2'])
    graph.add_edges(df_edge.loc[:, ['index_x','index_y']].values)
    if weight_col:
        graph.es["weight"] = df_edge.loc[:,weight_col].values
    if _extra_edge_col:
        graph.es[_extra_edge_col] = df_edge.loc[:,_extra_edge_col].values
    return graph

def Loops_Return_two_bins_no_dup(df_hic):
    ## Associated by promoter
    second_bin_columns = [3,4,5,0,1,2]+list(range(6,len(df_hic.columns),1))
    df_hic=df_hic.append(pd.DataFrame(df_hic.iloc[:, second_bin_columns].values, columns=df_hic.columns),sort=False).sort_index()
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

def graph_community_multilevel_Blondel(input_graph, cutoff):
    ## input graph should have at least one attribute: name
    df_vs = convert_graph_vs_to_df(input_graph)
    _col_vs_name='name'
    if (input_graph.is_weighted()):
        print ("Weighted Graph Cluster")
        structure = input_graph.community_multilevel(weights=input_graph.es['weight'] ,return_levels=False)
    else:
        structure = input_graph.community_multilevel(return_levels=False)
    df_vs['membership'] = structure.membership
    df_vs_cluster_group = df_vs.groupby('membership')
    
    ## Rank each cluster by number of bins
    cluster_name=[]
    cluster_num_vertices=[]
    for df_vs_cluster in df_vs_cluster_group:
        df_vs_inside_cluster = Cluster_Filter_by_Denisty(df_vs_cluster[1], _col_vs_name, 'degree', cutoff)
        #df_vs_inside_cluster =df_vs_cluster[1]
        df_cluster_coordiante = df_vs_inside_cluster[_col_vs_name].str.split(r"\:|-",expand=True)
        cluster_coordinate = 'chr'+df_cluster_coordiante.iloc[0,0]+':'+str(df_cluster_coordiante.iloc[:,1].astype(int).min())+'-'+str(df_cluster_coordiante.iloc[:,2].astype(int).max())
        cluster_name.append(cluster_coordinate) ##0: cluster name
        cluster_num_vertices.append(len(df_vs_inside_cluster)) # 1: num_vertices
    
    df_cluster_output = pd.DataFrame(data={'hub_name':cluster_name,'Num_vertices':cluster_num_vertices}).sort_values('Num_vertices', ascending=False)
    return df_cluster_output, df_vs_cluster_group

def Graph_Pagerank(_input_graph):
    input_graph = _input_graph
    input_graph.vs['pagerank'] = input_graph.pagerank(weights=input_graph.es['weight'])
    return input_graph

def Cluster_Filter_by_Denisty(_df_vs_cluster, _col_name, _core_col, _cutoff):
    ## Linear Denisty Threshold, 1 edge at least 1 anchor
    cutoff=_cutoff#0.5
    df_tem = _df_vs_cluster
    col_name='name'
    _core_col='pagerank'
    resolution=10000
    df_tem[col_name].str.split(r"\:|-",expand=True)
    df_tem = pd.concat( [df_tem[col_name].str.split(r"\:|-",expand=True),df_tem], axis=1)

    ## Define highest degree as summit
    num_core = 1
    core = df_tem.nlargest(int(num_core), _core_col).iloc[:,1].astype(int).mean()
    df_tem['density'] = df_tem['degree'].astype(float)/(abs(df_tem.iloc[:,1].astype(float)-float(core))**2)*resolution**2

    return df_tem[df_tem['density']>cutoff]

def graph_community_multilevel_Blondel_diff_level(input_graph, cutoff):
    ## input graph should have at least one attribute: name
    df_vs = convert_graph_vs_to_df(input_graph)
    _col_vs_name='name'
    
    if (input_graph.is_weighted()):
        print ("Weighted Graph Cluster")
        structure = input_graph.community_multilevel(weights=input_graph.es['weight'], return_levels=True)
    else:
        structure = input_graph.community_multilevel(return_levels=True)
    
    for tem_level in structure:
        print (tem_level.summary())
    df_vs['membership'] = structure[0].membership
    df_vs_cluster_group = df_vs.groupby('membership')
    
    ## Rank each cluster by number of bins
    cluster_name=[]
    cluster_num_vertices=[]
    for df_vs_cluster in df_vs_cluster_group:
        df_vs_inside_cluster = Cluster_Filter_by_Denisty(df_vs_cluster[1], _col_vs_name, 'degree', cutoff)
        if (len(df_vs_inside_cluster)>0):
            df_cluster_coordiante = df_vs_inside_cluster[_col_vs_name].str.split(r"\:|-",expand=True)
            #print (df_cluster_coordiante)
            cluster_coordinate = 'chr'+df_cluster_coordiante.iloc[0,0]+':'+str(df_cluster_coordiante.iloc[:,1].astype(int).min())+'-'+str(df_cluster_coordiante.iloc[:,2].astype(int).max())
            cluster_name.append(cluster_coordinate) ##0: cluster name
            cluster_num_vertices.append(len(df_vs_inside_cluster)) # 1: num_vertices
    
    df_cluster_output = pd.DataFrame(data={'hub_name':cluster_name,'Num_vertices':cluster_num_vertices}).sort_values('Num_vertices', ascending=False)
    return df_cluster_output, df_vs_cluster_group

def graph_community_multilevel_Blondel_diff_level_promoter(input_graph, cutoff):
    ## input graph should have at least one attribute: name
    df_vs = convert_graph_vs_to_df(input_graph)
    _col_vs_name='name'
    
    if (input_graph.is_weighted()):
        print ("Weighted Graph Cluster")
        structure = input_graph.community_multilevel(weights=input_graph.es['weight'], return_levels=True)
    else:
        structure = input_graph.community_multilevel(return_levels=True)
    
    for tem_level in structure:
        print (tem_level.summary())
    df_vs['membership'] = structure[0].membership
    df_vs_cluster_group = df_vs.groupby('membership')
    
    ## Rank each cluster by number of bins
    cluster_summary = []
    for df_vs_cluster in df_vs_cluster_group:
        df_cluster = df_vs_cluster[1] 
        if( len(df_cluster[df_cluster['Promoter']!=0])>0):
            for promoter_id in df_cluster[df_cluster['Promoter']!=0]['Promoter_gene_id'].unique():#[0]
                #promoter_id = df_cluster[df_cluster['Promoter']!=0]['Promoter_gene_id'].unique()[0]
                #print(promoter_id)
                if (promoter_id=='Myb'):
                    df_test_out = df_cluster
                df_vs_inside_cluster, cluster_coordinate = Cluster_Filter_by_Denisty_Promoter(df_cluster, _col_vs_name, promoter_id, cutoff)            
                cluster_summary.append( [cluster_coordinate, len(df_vs_inside_cluster), promoter_id])
                
    
    df_cluster_output = pd.DataFrame(data=cluster_summary, columns=['hub_name','Num_vertices','Promoter']).sort_values('Num_vertices', ascending=False)
    return df_cluster_output, df_vs_cluster_group, df_test_out

def Cluster_Filter_by_Denisty_Promoter(_df_vs_cluster, _col_name, _promoter_id, _cutoff):
    ## Linear Denisty Threshold, 1 edge at least 1 anchor
    cutoff=_cutoff#0.5
    df_tem = _df_vs_cluster
    col_name='name'
    promoter_id = _promoter_id
    resolution=10000
    df_tem[col_name].str.split(r"\:|-",expand=True)
    df_tem = pd.concat( [df_tem[col_name].str.split(r"\:|-",expand=True),df_tem], axis=1)
    
    ## Define Target Promoter as core
    core = df_tem[df_tem['Promoter_gene_id']==promoter_id].iloc[:,1:3].astype(int).sum(axis=1)/2
    #print (core)

    df_tem['density'] = df_tem['degree'].astype(float)/(abs(df_tem.iloc[:,1].astype(float)-float(core))**2)*resolution**2
    df_filtered_cluster_elements = df_tem[df_tem['density']>cutoff]
    
    df_cluster_coordiante = df_filtered_cluster_elements[col_name].str.split(r"\:|-",expand=True)
    
    cluster_coordinate = 'chr'+df_cluster_coordiante.iloc[0,0]+':'+str(df_cluster_coordiante.iloc[:,1].astype(int).min())+'-'+str(df_cluster_coordiante.iloc[:,2].astype(int).max())
    
    
    return df_filtered_cluster_elements, cluster_coordinate

def display_graph_vertex(input_graph, vs_idx_set):
    # Input graph and Input vertex index set
    for vs in graph_processed.vs.select(vs_idx_set):
        print ( "vs_idx:"+ str(vs.index)+ ' '+ str(vs.attributes()))
    return None

def display_graph_edge(input_graph, vs_idx_set):
    # Input graph and Input vertex index set
    for vs_idx in vs_idx_set:
        edges_from_vs = input_graph.es[input_graph.incident(vs_idx)]
        for es in edges_from_vs:
            print ( "es_idx:"+ str(es.index)+ ' '+ str(es.tuple))
    return None

def annotate_graph_with_feature_values_new(_input_graph, graph_name_col2bed, path_feature, feature_name, _feature_score, norm_factor=1.0):
    input_graph = _input_graph
    name_col2bed = graph_name_col2bed ## Default "name"
    Vs_Attrs_Name = feature_name ## such as 'Tcf1'
    if ( Vs_Attrs_Name not in input_graph.vs.attributes()):
        ## Convert vs to bed format in order to annotate
        df_vs_bed = convert_vs2bed(input_graph, name_col2bed)
        ### df_vs_bed to be annoted
        df_vs_bed.iloc[:,0]='chr'+ df_vs_bed.iloc[:,0].astype(str)
        Feature_vs = BedTool.from_dataframe(df_vs_bed).sort()

        PATH_Feature_A = path_feature ##
        df_A = pd.read_csv(PATH_Feature_A, sep="\t")
        Feature_A = BedTool.from_dataframe(df_A).sort()

        ## annotate A in vs
        Feature_vs_with_A = Feature_vs.intersect(Feature_A, wb=True, F=0.3) ## 30% maybe enough

        if (len(Feature_vs_with_A)>0):
            df_vs_with_A=pd.read_csv(Feature_vs_with_A.fn, sep="\t", names=df_vs_bed.columns.append(df_A.columns).values, header=None)
        else:
            df_vs_with_A=pd.DataFrame(columns=df_vs_bed.columns.append(df_A.columns))
        
        
        vs_score = _feature_score  ## 'such as logFC'
        vs_attrs_score = feature_name+'_'+vs_score
        input_graph.vs[Vs_Attrs_Name]=0
        input_graph.vs[vs_attrs_score]=0

        for df_vs in df_vs_with_A.groupby(name_col2bed): ### Default Define vertex attribute "name"
            input_graph.vs.select(name=df_vs[0])[Vs_Attrs_Name] = df_vs[1].shape[0]
            ### max Tcf1 binding
            if ( type(df_vs[1].loc[:,vs_score].head(1).values[0]) == str):
                input_graph.vs.select(name=df_vs[0])[vs_attrs_score] = df_vs[1].loc[:,vs_score].max()
            else:
                #print(df_vs[1].loc[:,vs_score])
                if (df_vs[1].shape[0]==1):
                    input_graph.vs.select(name=df_vs[0])[vs_attrs_score] = df_vs[1].loc[:,vs_score].values[0]/norm_factor
                else:
                    List_Feature = list(df_vs[1].loc[:,vs_score].values)
                    input_graph.vs.select(name=df_vs[0])[vs_attrs_score] = List2Str(List_Feature, norm_factor)#multiple save as str

        print ("Annotate " + Vs_Attrs_Name + " is finished.")
    else: 
        print ("Feature of " + Vs_Attrs_Name + " is already annoated. Skip.")
    
    return input_graph

def annotate_graph_with_feature_values(_input_graph, graph_name_col2bed, path_feature, feature_name, _feature_score):
    input_graph = _input_graph
    name_col2bed = graph_name_col2bed ## Default "name"
    Vs_Attrs_Name = feature_name ## such as 'Tcf1'
    if (Vs_Attrs_Name not in _input_graph.vs.attributes()):
        ## Convert vs to bed format in order to annotate
        df_vs_bed = convert_vs2bed(input_graph, name_col2bed)
        ### df_vs_bed to be annoted
        Feature_vs = BedTool.from_dataframe(df_vs_bed).sort()

        PATH_Feature_A = path_feature ##
        df_A = pd.read_csv(PATH_Feature_A, sep="\t")
        Feature_A = BedTool.from_dataframe(df_A).sort()

        ## annotate A in vs
        Feature_vs_with_A = Feature_vs.intersect(Feature_A, wb=True, F=0.3)

        if (len(Feature_vs_with_A)>0):
            df_vs_with_A=pd.read_csv(Feature_vs_with_A.fn, sep="\t", names=df_vs_bed.columns.append(df_A.columns).values, header=None)
        else:
            df_vs_with_A=pd.DataFrame(columns=df_vs_bed.columns.append(df_A.columns))
        
        
        vs_score = _feature_score  ## 'such as logFC'
        vs_attrs_score = Vs_Attrs_Name+'_'+vs_score
        input_graph.vs[Vs_Attrs_Name]=0
        input_graph.vs[vs_attrs_score]=0
        for df_vs in df_vs_with_A.groupby(name_col2bed): ### Default Define vertex attribute "name"
            input_graph.vs.select(name=df_vs[0])[Vs_Attrs_Name] = df_vs[1].shape[0]
            ### max Tcf1 binding
            if ( type(df_vs[1].loc[:,vs_score].head(1).values[0]) == str):
                input_graph.vs.select(name=df_vs[0])[vs_attrs_score] = df_vs[1].loc[:,vs_score].max()
            else:
                #print(df_vs[1].loc[:,vs_score])
                input_graph.vs.select(name=df_vs[0])[vs_attrs_score] = df_vs[1].loc[:,vs_score].mean()
        print ("Annotate " + Vs_Attrs_Name + " is finished.")
    else: 
        print ("Feature of " + Vs_Attrs_Name + " is already annoated. Skip.")
    
    return input_graph

def Return_Graph_of_Gene(_input_graph, _gene, _search_depth):
    search_depth=_search_depth
    graph_input = _input_graph
    gene_request = _gene
    vertex_set=set()
    if ( len(graph_input.vs.select(Promoter_gene_id=gene_request))>0 ):
        vertex_set.add(graph_input.vs.select(Promoter_gene_id=gene_request)[0].index)
        final_subgrapph=None
        for i in range(search_depth):
            for vertex_index in list(vertex_set):
                graph_select_edges = graph_input.es[graph_input.incident(vertex_index)]
                for edge_in_graph in graph_select_edges:
                    vertex_set.add(edge_in_graph.tuple[0])
                    vertex_set.add(edge_in_graph.tuple[1])

            #print (vertex_set)
            final_subgrapph = graph_input.induced_subgraph(vertex_set)
    else:
        print ('Gene Not Included!')
        final_subgrapph=None
    
    return final_subgrapph

def Read_Interaction(_PATH_interaction, _resolution, _col_fore, _col_back):
    PATH_interaction=_PATH_interaction
    col_fore = _col_fore
    col_back  = _col_back
    resolution = _resolution
    
    df_interaction = pd.read_csv(PATH_interaction, sep="\t").fillna(0)
    df_interaction = df_interaction[df_interaction.iloc[:,1]!=df_interaction.iloc[:,2]] ### remove self interaction
    df_interaction.loc[:,'#chr']=df_interaction.iloc[:,0].replace('chr','')
    df_interaction.loc[:,'#chr1']=df_interaction.iloc[:,0]
    df_interaction.loc[:,'x1']=df_interaction.iloc[:,1].astype(int)
    df_interaction.loc[:,'x2']=df_interaction.iloc[:,1].astype(int)+int(resolution)
    df_interaction.loc[:,'chr2']=df_interaction.iloc[:,0]
    df_interaction.loc[:,'y1']=df_interaction.iloc[:,2].astype(int)
    df_interaction.loc[:,'y2']=df_interaction.iloc[:,2].astype(int)+int(resolution)

    df_interaction.loc[:,'log_FC'] = np.log2(df_interaction.loc[:,col_fore].replace(0,0.1) / df_interaction.loc[:,col_back].replace(0,0.1) )
    #df_interaction.loc[:,'GeneID'] = "id_"+df_interaction.index.astype(str)
    df_interaction = df_interaction.loc[:,['#chr1','x1','x2','chr2','y1','y2','log_FC', col_fore, col_back]]
    return df_interaction

def calculate_pvalue_for_hub(_PATH_interaction, _df_Hubs, _col_fore, _col_back):
    ## Calculate pvalue for each hub
    PATH_interaction = _PATH_interaction
    col_fore = _col_fore
    col_back  = _col_back
    df_Hub_top = convert_cluster2bed(_df_Hubs, 'hub_name').reset_index().drop('index', axis=1)
    
    ## Associated each Hub with interaction and pvalue
    ########################################################################################################
    df_inter = pd.read_csv(PATH_interaction, sep="\t").fillna(0)
    df_inter = df_inter[df_inter.iloc[:, 1]!=df_inter.iloc[:, 2]]
    df_inter.loc[:,'#chr']= 'chr'+df_inter.iloc[:,0].astype(str)
    Feature_interaction = BedTool.from_dataframe(df_inter).sort()
    Feature_hub = BedTool.from_dataframe(df_Hub_top).sort()
    ########################################################################################################
    ## calculate all interactions inside a hub
    Feature_Hub_interaction = Feature_hub.intersect(Feature_interaction, wa=True, wb=True, F=1.0)
    col_name = df_Hub_top.columns.append(df_inter.columns)
    df_Feature_Hub_interaction = pd.read_csv(Feature_Hub_interaction.fn, sep='\t', names=col_name)
    df_Feature_Hub_interaction_group = df_Feature_Hub_interaction.groupby('hub_name')

    ########################################################################################################
    ### calculate a pvalue for each hub
    hub_sum=[]
    for hub in df_Feature_Hub_interaction_group:
        #print (hub[0])
        df_hub = hub[1]
        data_for_test = df_hub.loc[:, col_back]  - df_hub.loc[:, col_fore]
        w, pvalue_hub = stats.wilcoxon(data_for_test)#, alternative='less')
        hub_sum.append([hub[0], df_hub.Num_vertices.unique()[0], pvalue_hub])
        #break

    df_hub_summary = pd.DataFrame( data = hub_sum, columns=['hub_name', 'Num_vertices', 'pvalue'])
    df_hub_summary = df_Hub_top.merge(df_hub_summary, on=['hub_name','Num_vertices'], how='inner').sort_values(by='pvalue')
    
    return df_hub_summary
### End of Visulization
####################################################################################
### FUNCTION
### FUNCTIONS
def main(argv):
	desc="Collect HiC Interaction in txt format, rank interaction change Hub. Input Format should be: #chr	bin1	bin2	Cond1	Cond2"
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Path to Input HiC file in txt format", metavar="<file>")
	parser.add_option("-f", "--foreground_name", action="store", type="string",
			dest="fore_name", help="Name of condition as foreground.", metavar="<str>")
	parser.add_option("-b", "--background_name", action="store", type="string",
			dest="back_name", help="Name of condition as background.", metavar="<str>")
	parser.add_option("-r", "--resolution", action="store", type="int",
		dest="res", help="Resolution of HiC txt", metavar="<int>")
	parser.add_option("-d", "--filtered_density", action="store", type="float",
		dest="density_cutoff", help="Density cutoff for hub shriking.", metavar="<float>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 5:
		parser.print_help()
		sys.exit(1)
	
	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in txt format: %s" % opt.input_path)
	print ("Foreground Condition: %s" % opt.fore_name)
	print ("Background Condition: %s" % opt.back_name)
	print ("Resolution %i" % opt.res)
	print ("Density cutoff for hub shriking: %f" % opt.density_cutoff)
	print ("End of Summary.")
	print (" ")
	
	## parameters
	PATH_INPUT=opt.input_path
	col_fore = opt.fore_name
	col_back  = opt.back_name
	resolution = opt.res
	density_cutoff = opt.density_cutoff


#### Main 
	df_hic = Read_Interaction(PATH_INPUT, resolution, col_fore, col_back)

	for diff_hub in [col_fore, col_back]:
		if (diff_hub==col_fore):
			df_diff_hic  = df_hic[(df_hic['log_FC']>0)]
		else:
			df_diff_hic  = df_hic[(df_hic['log_FC']<0)]
		print ("Processing Data For Diff Hub On: "+diff_hub)
		g_graph = convert_loops_to_graph(df_diff_hic, col_back, col_fore)
		graph_processed = Graph_Pagerank(g_graph)
		df_Diff_Hub, df_Hub_Groups = graph_community_multilevel_Blondel_diff_level(graph_processed, density_cutoff)
		df_hub_summary = calculate_pvalue_for_hub(PATH_INPUT, df_Diff_Hub, col_fore, col_back)
		df_hub_summary.to_csv(str(len(df_hub_summary))+'_'+diff_hub+'_Diff_hub.txt',sep='\t', index=None)
		print("Result is saved as: " + str(len(df_hub_summary))+'_'+diff_hub+'_Diff_hub.txt')
		print(" ")
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
