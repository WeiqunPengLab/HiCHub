import pandas as pd
import numpy as np
import igraph as ig
#import pybedtools
from scipy import stats
from optparse import OptionParser
import sys, os, multiprocessing
import random

random.seed(10)

def revise_hub(hubs):
    df_test = hubs
    #df_test = df_test[df_test['-log10(pvalue)']>20]
    
    #pyramid = pd.DataFrame()
    x = df_test[df_test['left_hub_anchor'] == df_test['right_hub_anchor']]
    pyramid = pd.DataFrame(data={'0':x.loc[:,'left_hub_anchor'].str.split(':', expand = True)[0],
                                 '1':x.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int),
                                 '2':x.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int),
                                 'left_hub_anchor':x.loc[:,'left_hub_anchor'],
                                 'right_hub_anchor':x.loc[:,'right_hub_anchor'],
                                 '-log10(pvalue)':x.loc[:,'-log10(pvalue)']})

    y = df_test[df_test['left_hub_anchor'] != df_test['right_hub_anchor']]
    stripe_1 = pd.DataFrame(data={'0':y.loc[:,'left_hub_anchor'].str.split(':', expand = True)[0],
                                  '1':y.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int),
                                  '2':y.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int),
                                  'left_hub_anchor':y.loc[:,'left_hub_anchor'],
                                  'right_hub_anchor':y.loc[:,'right_hub_anchor'],
                                  '-log10(pvalue)':y.loc[:,'-log10(pvalue)']})
    #stripe_1['0'] = stripe_1.loc[:,'left_hub_anchor'].str.split(':', expand = True)[0]
    #stripe_1['1'] = stripe_1.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    #stripe_1['2'] = stripe_1.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    #stripe_2 = pd.DataFrame()
    z = df_test[df_test['left_hub_anchor'] != df_test['right_hub_anchor']]
    stripe_2 = pd.DataFrame(data={'0':z.loc[:,'left_hub_anchor'].str.split(':', expand = True)[0],
                                  '1':z.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int),
                                  '2':z.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int),
                                  'left_hub_anchor':z.loc[:,'left_hub_anchor'],
                                  'right_hub_anchor':z.loc[:,'right_hub_anchor'],
                                  '-log10(pvalue)':z.loc[:,'-log10(pvalue)']})
    #stripe_2['0'] = stripe_2.loc[:,'right_hub_anchor'].str.split(':', expand = True)[0]
    #stripe_2['1'] = stripe_2.loc[:,'right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    #stripe_2['2'] = stripe_2.loc[:,'right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)
    
    '''
    pyramid.loc[:,'0'] = pyramid.loc[:,'left_hub_anchor'].str.split(':', expand = True)[0]
    pyramid.loc[:,'1'] = pyramid.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    pyramid.loc[:,'2'] = pyramid.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    #stripe_1 = pd.DataFrame()
    stripe_1 = df_test[df_test['left_hub_anchor'] != df_test['right_hub_anchor']]
    stripe_1.loc[:,'0'] = stripe_1.loc[:,'left_hub_anchor'].str.split(':', expand = True)[0]
    stripe_1.loc[:,'1'] = stripe_1.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    stripe_1.loc[:,'2'] = stripe_1.loc[:,'left_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    #stripe_2 = pd.DataFrame()
    stripe_2 = df_test[df_test['left_hub_anchor'] != df_test['right_hub_anchor']]
    stripe_2.loc[:,'0'] = stripe_2.loc[:,'right_hub_anchor'].str.split(':', expand = True)[0]
    stripe_2.loc[:,'1'] = stripe_2.loc[:,'right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    stripe_2.loc[:,'2'] = stripe_2.loc[:,'right_hub_anchor'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)
    '''
    
    stripe = pd.DataFrame()
    stripe = pd.concat([stripe, stripe_1], axis=0)#stripe.append(stripe_1)
    stripe = pd.concat([stripe, stripe_2], axis=0)#stripe = stripe.append(stripe_2)
    stripe.drop_duplicates(subset=['0','1','2'],keep='first', inplace=True)


    total = pd.DataFrame()
    total = pd.concat([total, stripe], axis=0)#total = total.append(stripe)
    total = pd.concat([total, pyramid], axis=0)#total = total.append(pyramid)
    total.drop_duplicates(subset=['0','1','2'],keep='first', inplace=True)
    total['label'] = 'H1ESC'

    total = total.loc[:,['0','1','2','left_hub_anchor','right_hub_anchor','-log10(pvalue)','label']]
    #total.to_csv('top_H1ESC_regions.bed', sep = '\t', index = None)
    return total

def find_hub(gene_name, input_path, file_name, file_label):
    
    promoter = pd.read_csv('promoter.bed', sep='\t')
    
    promoter = promoter[promoter['gene_id'] == gene_name]
    if len(promoter) == 1:
        chrom = promoter.iloc[0,0]
        start = promoter.iloc[0,1]
    else:
        print('Sorry, no such gene, please enter the right gene name.')
        return 'x'
    
    k=0
    for i in range(len(file_name)):
        hub = pd.read_csv(file_name[i], sep='\t')
        hub = revise_hub(hub)
        hub = hub[hub['0'] == chrom]
        hub = hub[hub['1'] <= start]
        hub = hub[hub['2'] >= start]
        hub = hub.loc[:,['left_hub_anchor','right_hub_anchor','-log10(pvalue)']]
        if len(hub) > 0:
            print('The gene ' + gene_name + ' is in '+file_label[i]+' hubs:')
            ttt = i
            k = k+1
            x=[]
            for i in range(len(hub)):
                x.append(i+1)
            print(hub.to_string(index=False))
    if k == 0:
        print('The '+gene_name+' does not exist in any hubs.')
        return 'x'
    return ttt
	
	

def Add_pesudo_count(_df_test, _count, _col1, _col2):
    df_test = _df_test.copy()
    count = _count
    col1 = _col1
    col2 = _col2
    #print (df_test)
        
    df_test[col1] = df_test[col1] + count
    df_test[col2] = df_test[col2] + count
    #print (df_test)
    return df_test

def LOESS_Norm_df (df_test, _col1, _col2):
    ## this is a similar approach as LOESS Normalization
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
        df_out = pd.concat([df_out, df_bbb], axis=0)
        #df_out = df_out.append(df_bbb)
        df_out = df_out.sort_index()
    return df_out.sort_index()

def plot_igragh(edges, cluster, gene):
    test = cluster#pd.read_csv('cluster_annotated.bed',sep='\t')#,names=['chr','node','cluster'])
    #test = test[test['gene_id'] == 'LSR']
    test = test.replace('up','triangle-up')
    test = test.replace('down','triangle-down')
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
            shape_dict.append('circle')
        elif test.iloc[i,4] == gene:
            shape_dict.append('circle')
        elif test.iloc[i,4] =='0' and test.iloc[i,5] != '0':
            shape_dict.append(test.iloc[i,5])
        elif test.iloc[i,4] !='0' and test.iloc[i,5] != '0':
            shape_dict.append('circle')
        else:
            shape_dict.append('circle')
    visual_style["vertex_shape"] = shape_dict
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
    visual_style['edge_width'] = [0.3]
    visual_style['label_dist'] = 20

    # 图尺寸
    #visual_style["bbox"] = (600, 480)
    # 边缘距离
    #visual_style["margin"] = 50
    # 布局方式
    #visual_style["layout"] = layout
    # -----------------------画图-----------------------------
    ig.plot(graph_term, **visual_style)

    #ig.plot(g, layout = 'kk').show()
    return None



def find_cluster(input_path, cluster_path, gene_name, fore_name, back_name, FC, cut_off, promoter_path):
	
    test = pd.read_csv(input_path, sep='\t', dtype={'#chr':str}, iterator=True, chunksize=100000)


    out = pd.DataFrame()
    for chunk in test:
        chunk = chunk.fillna(0)
        chunk = chunk[chunk[fore_name]+chunk[back_name] > cut_off]
        out = pd.concat([out,chunk],axis=0)
			
    edges_file = out
	
    cluster_file = pd.read_csv(cluster_path, sep='\t')
    promoter = pd.read_csv(promoter_path, sep='\t')
    
    promoter = promoter[promoter['gene_id'] == gene_name]
    if len(promoter) == 1:
        chrom = promoter.iloc[0,0]
        start = promoter.iloc[0,1]
    else:
        print('Sorry, no such gene, please enter the right gene name.')
        sys.exit(1)
    cluster = cluster_file[cluster_file['#chr'] == chrom]
    cluster = cluster[cluster['start'] <= start]
    cluster = cluster[cluster['end'] >= start]
    
    if len(cluster) == 1:    
        cluster_num = cluster.iloc[0,3]
    else:
        print('Sorry, this gene is not in any of our clusters.')
        sys.exit(1)
    cluster_input = cluster_file[cluster_file['#chr'] == chrom]
    cluster_input = cluster_input[cluster_input['cluster'] == cluster_num]
    
    edges_file = edges_file['chr'+edges_file['#chr'] == chrom]
    edges_input = LOESS_Norm_df( edges_file, fore_name, back_name)
    edges_input['diff'] = edges_input[fore_name] / edges_input[back_name]
    edges_input = edges_input[edges_input['diff']>FC]
    edges_input = edges_input.loc[:,['#chr','bin1','bin2','diff']]
    
    plot_igragh(edges_input, cluster_input, gene_name)
    return None
    
    

####################################################################################
### FUNCTION
### FUNCTIONS
def run(opt):
	## parameters
	## parameters
	PATH_INPUT=opt.input_path
	col_fore = opt.label.split(',')[0]
	col_back = opt.label.split(',')[1]
	cut_off = opt.cut_off_value
	FC = opt.foldchange
	gene  = opt.gene_name.split(',')
	promoter_path = opt.promoter_path
	

	print (" ")
	print("Run main")
	print ("Here is the Summary of your input.")	
	print ("Input Path of HiC file in txt format: %s" % PATH_INPUT)
	print ("FC for qualified difference is: %s" % FC)
	print ("Threshold value for sum of row counts is: %s" % cut_off)
	print ("First label is: %s" % col_fore)
	print ("Second label is: %s" % col_back)	
	for i in gene:
		print ("Name of gene: %s" % i)
	print ("Path of promoter: %s" % promoter_path)
	print ("End of Summary.")
	print (" ")
	print (" ")


#### Main
	for sub_name in gene:
		File_Name_Sets = [col_fore+'_specific_regions.bed', col_back+'_specific_regions.bed']
		File_Label_Sets = [col_fore,col_back]
		#print(File_Label_Sets)
		result1 = find_hub(sub_name, 'no_need', File_Name_Sets, File_Label_Sets)

		if result1 != 'x':
			cluster = 'cluster_annotated_'+ File_Label_Sets[result1] + '.txt'
			if result1 == 0:
				find_cluster(PATH_INPUT, cluster, sub_name, File_Label_Sets[0], File_Label_Sets[1], FC, cut_off, promoter_path)
			if result1 == 1:
				find_cluster(PATH_INPUT, cluster, sub_name, File_Label_Sets[1], File_Label_Sets[0], FC, cut_off, promoter_path)
			
	return None

def main(argv):
	
	desc="Please enter the file of hubs, clusters and original .txt file."
	parser = OptionParser(description=desc)
	
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Input .txt file you want to call hubs..", metavar="<file>")
	parser.add_option("-d", "--FC", action="store", type="float", default = 1.0,
		dest="foldchange", help="Optional: FC to evaluate qualified difference, default=1.0.", metavar="<float>")
	parser.add_option("-c", "--cut_off", action="store", type="float", default = 10,
		dest="cut_off_value", help="Optional: Remove sum of row counts smaller than this value, default=10.", metavar="<float>")
	
	parser.add_option("-l", "--label", action="store", type="string",
			dest="label", help="Name of labels you have used to call hubs. Delimiter ',', for example 'label1,label2'.", metavar="<str>")
	
	parser.add_option("-p", "--promoter", action="store", type="string",
			dest="promoter_path", help="File name of gene promoters.", metavar="<file>")
	
	parser.add_option("-n", "--gene name", action="store", type="string",
			dest="gene_name", help="Name of gene as input.", metavar="<str>")


    
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	## parameters
	PATH_INPUT=opt.input_path
	col_fore = opt.label.split(',')[0]
	col_back = opt.label.split(',')[1]
	cut_off = opt.cut_off_value
	FC = opt.foldchange
    
	gene  = opt.gene_name.split(',')
    
	promoter_path = opt.promoter_path
	

	print (" ")
	print("Run main")
	print ("Here is the Summary of your input.")	
	print ("Input Path of HiC file in txt format: %s" % PATH_INPUT)
	print ("FC for qualified difference is: %s" % FC)
	print ("Threshold value for sum of row counts is: %s" % cut_off)
	print ("First label is: %s" % col_fore)
	print ("Second label is: %s" % col_back)
	for i in gene:
		print ("Name of gene: %s" % i)
	print ("Path of promoter: %s" % promoter_path)
	print ("End of Summary.")
	print (" ")
	print (" ")
	


#### Main
	for sub_name in gene:
		File_Name_Sets = [col_fore+'_specific_regions.bed', col_back+'_specific_regions.bed']
		File_Label_Sets = [col_fore,col_back]
		#print(File_Label_Sets)
		result1 = find_hub(sub_name, 'no_need', File_Name_Sets, File_Label_Sets)

		if result1 != 'x':
			cluster = 'cluster_annotated_'+ File_Label_Sets[result1] + '.txt'
			if result1 == 0:
				find_cluster(PATH_INPUT, cluster, sub_name, File_Label_Sets[0], File_Label_Sets[1], FC, cut_off, promoter_path)
			if result1 == 1:
				find_cluster(PATH_INPUT, cluster, sub_name, File_Label_Sets[1], File_Label_Sets[0], FC, cut_off, promoter_path)


#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
