import pandas as pd
import numpy as np
import igraph as ig
#import pybedtools
from scipy import stats
from optparse import OptionParser
import sys, os, multiprocessing
import random
import matplotlib.pyplot as plt

random.seed(5)


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
       
        df_out = pd.concat([df_out,df_bbb], axis=0)#df_out = df_out.append(df_bbb)
        df_out = df_out.sort_index()
    return df_out.sort_index()

def plot_igragh(edges, cluster, gene):
    test = cluster#pd.read_csv('cluster_final.bed',sep='\t')#,names=['chr','node','cluster'])
    #test = test[test['gene_id'] == 'LSR']
    test = test.replace('up',r"^")
    test = test.replace('down',r"v")
    #test = test[test['signal'] == 'triangle-down']
    test['promoter'] = test['gene_id'].apply(lambda x: 0 if x == '0' else(30))
    #test = test[test['start']==20790000]
    #test = test[test['cluster']=='cluster330']
    #test = test.replace('0',np.nan)
    df_bins = test.loc[:,['start']].reset_index()
    df_bins = df_bins.loc[:,['start']]
    df_bins['label'] = 'N'
    for i in range(len(df_bins)):
        if df_bins.loc[i,'start'] == 41190000:
            df_bins.loc[i,'label'] = 'promoter'
            

    Num_vs = len(df_bins)
    graph_term = ig.Graph()

    graph_term.add_vertices(df_bins['start'].astype(str))
    graph_term.vs['l'] = df_bins.loc[:,'label']
    #graph_term.vs['color'] = test.loc[:,'color']
    #graph_term.vs['name'] = df_bins.loc[:,'node'].astype(str)

    df_edge = edges#pd.read_csv('chr19_H1ESC.txt', sep='\t', dtype={'#chr':str})
    df_edge = df_edge[df_edge['bin1'].isin(df_bins['start'])]
    df_edge = df_edge[df_edge['bin2'].isin(df_bins['start'])]
    edges = []
    #n = df_edge.loc[:,['bin1','bin2']].values

    for i in range(len(df_edge)):
        edges.append((df_edge.iloc[i,1].astype(str),df_edge.iloc[i,2].astype(str)))

    graph_term.add_edges(edges)
    #print(graph_term.summary())

    #structure = graph_term.community_multilevel(return_levels=True)
    #graph_term = structure[0].subgraph(5)
    #print(graph_term.vs['name'])

    #test.to_csv('plot.bed',sep='\t',index=None)
    # -----------------------设置参数-----------------------------

    # 参数集合。visual_style是一个参数字典，可以动态添加想要个性化设定的参数
    visual_style = {}
    # 根据相对面积，设置点的大小
    #size_dict = {'promoter':30,'N':0}
    size_dict = []
    for i in range(len(test)):
        if test.iloc[i,4] !='0' and test.iloc[i,5] == '0' and test.iloc[i,4] != gene:
            size_dict.append(25*1.5)
        elif test.iloc[i,4] == gene:
            size_dict.append(35*1.5)
        elif test.iloc[i,4] =='0' and test.iloc[i,5] != '0':
            size_dict.append(15*1.5)
        elif test.iloc[i,4] !='0' and test.iloc[i,5] != '0':
            size_dict.append(25*1.5)
        else:
            size_dict.append(0)
    visual_style["vertex_size"] = size_dict
    # 根据国家实力，设置点的颜色
    color_dict = []
    for i in range(len(test)):
        if test.iloc[i,4] !='0' and test.iloc[i,4] != gene and test.iloc[i,5] == '0':
            color_dict.append('orange')
        elif test.iloc[i,4] == gene:
            color_dict.append('orange')
        elif test.iloc[i,4] =='0' and test.iloc[i,5] == 'triangle-up':
            color_dict.append('green')
        elif test.iloc[i,4] =='0' and test.iloc[i,5] == 'triangle-down':
            color_dict.append('red')
        elif test.iloc[i,4] !='0' and test.iloc[i,5] != '0':
            color_dict.append('orange')
        else:
            color_dict.append('orange')
        #color_dict.append(test.iloc[i,7])
    visual_style["vertex_color"] = color_dict
    #shape_dict = {'promoter':'triangle-down','N':'circle'}

    shape_dict = []
    for i in range(len(test)):
        if test.iloc[i,4] !='0' and test.iloc[i,4] != gene and test.iloc[i,5] == '0':
            shape_dict.append(r".")
        elif test.iloc[i,4] == gene:
            shape_dict.append(r".")
        elif test.iloc[i,4] =='0' and test.iloc[i,5] != '0':
            shape_dict.append(test.iloc[i,5])
        elif test.iloc[i,4] !='0' and test.iloc[i,5] != '0':
            shape_dict.append(r".")
        else:
            shape_dict.append(r".")
    #visual_style["vertex_shape"] = ['.']
    #['triangle-up','triangle-up','triangle-up','triangle-up','triangle-up','triangle-down',
                                    #'triangle-down','triangle-up','triangle-down','triangle-down','triangle-up','triangle-down',
                                    #'triangle-down','triangle-down',
                                    #'triangle-up','triangle-down','triangle-up','triangle-down','triangle-up','triangle-down','triangle-up','triangle-down',
                                    #'triangle-up','triangle-up','triangle-up','triangle-up',]

    label_dict = []
    for i in range(len(test)):
        if test.iloc[i,4] !='0':    
            label_dict.append(test.iloc[i,4])
        else:
            label_dict.append(None)
    #label_dict = {'promoter':'Promoter','N': None}
    visual_style["vertex_label"] = label_dict
    visual_style["label_size"] = 20
    # 边的粗细（这里随机生成）
    visual_style['edge_width'] = 0.3
    visual_style['label_dist'] = 20

    # 图尺寸
    #visual_style["bbox"] = (600, 480)
    # 边缘距离
    #visual_style["margin"] = 50
    # 布局方式
    #visual_style["layout"] = layout
    # -----------------------画图-----------------------------
    fig, ax = plt.subplots()
    ig.plot(graph_term, target = ax, **visual_style)
    plt.savefig('result.png')
    plt.show()

    #ig.plot(g, layout = 'kk').show()
    return None

def find_cluster(input_path, cluster_path, chrom, node, fore_name, back_name, _cut_off, FC, gene_name):
	
	test = pd.read_csv(input_path, sep='\t', dtype={'#chr':str}, iterator=True, chunksize=100000)

	out = pd.DataFrame()
	for chunk in test:
		chunk = chunk.fillna(0)
		chunk = chunk[chunk[fore_name]+chunk[back_name] > _cut_off]
		out = pd.concat([out,chunk],axis=0)
        
	edges_file = out
	edges_file = edges_file['chr'+ edges_file['#chr'] == chrom]
	
	#edges_file = pd.read_csv(input_path, sep='\t', dtype={'#chr':str})
	cluster_file = pd.read_csv(cluster_path, sep='\t')
    #promoter = pd.read_csv('promoter.bed', sep='\t')
    
    #promoter = promoter[promoter['gene_id'] == gene_name]
    #if len(promoter) == 1:
	chrom = chrom
	start = node
	print(chrom, start)
    #else:
        #print('Sorry, no such gene, please enter the right gene name.')
        #sys.exit(1)
	cluster = cluster_file[cluster_file['#chr'] == chrom]
	cluster = cluster[cluster['start'] == start]
	#cluster = cluster[cluster['end'] >= start]
	print(cluster)
    
	if len(cluster) == 1:    
		cluster_num = cluster.iloc[0,3]
	else:
		print('Sorry, this node is not in any of our clusters.')
		sys.exit(1)
	cluster_input = cluster_file[cluster_file['#chr'] == chrom]
	cluster_input = cluster_input[cluster_input['cluster'] == cluster_num]
    
	edges_file = edges_file['chr'+edges_file['#chr'] == chrom]
	edges_input = LOESS_Norm_df( edges_file, fore_name, back_name)
	edges_input['FC'] = edges_input[fore_name]/edges_input[back_name]
	edges_input = edges_input[edges_input['FC']> FC]
	edges_input = edges_input.loc[:,['#chr','bin1','bin2','FC']]
    
	plot_igragh(edges_input, cluster_input, gene_name)
	return None

def main(argv):
	
	desc="Collect HiC Interaction in txt format, rank interaction change Hub. Input Format should be: #chr	bin1	bin2	Cond1	Cond2"
	parser = OptionParser(description=desc)
	
	
	
	parser.add_option("-p", "--chrom_position", action="store", type="string",
			dest="chr_position", help="The chrom position for your plot.", metavar="<file>")
	parser.add_option("-n", "--node_position", action="store", type="int",
			dest="node_position", help="The node position for your plot.", metavar="<file>")
	
	parser.add_option("-g", "--gene_name", action="store", type="string", default = 'NONE',
			dest="gene_name", help="The node position for your plot.", metavar="<file>")

	parser.add_option("-o", "--cluster_files", action="store", type="string",
			dest="cluster_path", help="Cluster information that you had.", metavar="<str>")
	
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Path to Input HiC file in txt format", metavar="<file>")
	parser.add_option("-f", "--foreground_name", action="store", type="string",
			dest="fore_name", help="Name of condition as foreground.", metavar="<str>")
	parser.add_option("-b", "--background_name", action="store", type="string",
			dest="back_name", help="Name of condition as background.", metavar="<str>")

	
	parser.add_option("-d", "--FC", action="store", type="float", default = 1.0,
		dest="foldchange", help="Optional: FC to evaluate qualified difference, default=1.0", metavar="<float>")
	parser.add_option("-c", "--cut_off", action="store", type="float", default = 10,
		dest="cut_off_value", help="cut-off value for sum of row counts, default=10", metavar="<float>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
		
	PATH_INPUT=opt.input_path
	cluster = opt.cluster_path
	fore_name = opt.fore_name
	back_name = opt.back_name
	cut_off = opt.cut_off_value
	FC = opt.foldchange
	chrom = opt.chr_position
	node = opt.node_position
	gene_name = opt.gene_name
	
	print (" ")
	print ("Run main")
	print ("Here is the Summary of your input.")
	
	print ("Input Path of HiC file in txt format: %s" % PATH_INPUT)
	print ("Input path of cluster file: %s" % cluster)
	print ("Input chrom: %s" % chrom)
	print ("Input node: %s" % node)
	print ("Enphisize gene name: %s" % gene_name)
	print ("FC of MA normaliztion: %s" % FC)
	print ("Cut off value of MA normaliztion: %s" % cut_off)
	print ("Foreground Condition: %s" % fore_name)
	print ("Background Condition: %s" % back_name)
	print ("End of Summary.")
	print (" ")
	
	find_cluster(PATH_INPUT, cluster, chrom, node, fore_name, back_name, cut_off, FC, gene_name)
	

def run(argv):
	## parameters
	
	PATH_INPUT=argv.input_path
	cluster = argv.cluster_path
	fore_name = argv.fore_name
	back_name = argv.back_name
	cut_off = argv.cut_off_value
	FC = argv.foldchange
	chrom = argv.chr_position
	node = argv.node_position
	gene_name = argv.gene_name
	
	print (" ")
	print ("Run main")
	print ("Here is the Summary of your input.")
	
	print ("Input Path of HiC file in txt format: %s" % PATH_INPUT)
	print ("Input path of cluster file: %s" % cluster)
	print ("Input chrom: %s" % chrom)
	print ("Input node: %s" % node)
	print ("Enphisize gene name: %s" % gene_name)
	print ("FC of MA normaliztion: %s" % FC)
	print ("Cut off value of MA normaliztion: %s" % cut_off)
	print ("Foreground Condition: %s" % fore_name)
	print ("Background Condition: %s" % back_name)
	print ("End of Summary.")
	print (" ")
	
#### Main 
	find_cluster(PATH_INPUT, cluster, chrom, node, fore_name, back_name, cut_off, FC, gene_name)
	
	print (" ")
	return None

if __name__ == "__main__":
	main(sys.argv)
	