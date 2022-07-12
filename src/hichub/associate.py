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

def revise1(attach_name, PATH):
    test = pd.read_csv('cluster.txt', sep='\t', names=['chr','node','cluster'], dtype={'chr':str,'start':int,'end':int})

    test['#chr'] = 'chr' + test['chr']
    test['start'] = test['node']
    test['end'] = test['node'] + 10000

    test = test.loc[:,['#chr','start','end','cluster']]

    test.to_csv('cluster_combine.bed', sep='\t', index=None)

    test1 = pd.read_csv(PATH + attach_name, sep='\t', dtype={'chr':str,'start':int,'end':int})

    test1['mid'] = (0.5*(test1['start'] + test1['end'])).astype(int)

    test1['start'] = test1['mid']
    test1['end'] = test1['mid']+1

    test1 = test1.loc[:,['#chr','start','end','signal']]

    test1.to_csv('ATAC_combine.bed', sep='\t', index=None)
    return None

def revise2(PATH, promoter_name):
    DNase = pybedtools.BedTool('ATAC_combine.bed')
    promoter = pybedtools.BedTool(promoter_name)
    cluster = pybedtools.BedTool('cluster_combine.bed')

    c1 = promoter.intersect(cluster, wa=True, wb=True)

    c1.saveas('P_C.bed')

    p = pd.read_csv('P_C.bed', sep='\t', names=['#chr_x','start_x','end_x','gene_id','#chr_y','start_y','end_y','cluster_y'])
    c = pd.read_csv('cluster_combine.bed', sep='\t')

    p = p.loc[:,['#chr_y','start_y','end_y','cluster_y','gene_id']]

    c = c.merge(p, left_on='start', right_on='start_y', how='outer')

    c = c.loc[:,['#chr','start','end','cluster','gene_id']]

    c.to_csv('cluster_with_promoter.bed', sep='\t', index=None)
    return None

def revise3():
    test = pd.read_csv('ATAC_combine.bed', sep='\t')

    test['signal'] = test['signal']#.apply(lambda x: 'up' if x >= 0 else('down'))
    test['start'] = test['start'].astype(int)
    test['end'] = test['end'].astype(int)
    test = test.loc[:,['#chr','start','end','signal']]

    test.to_csv('ATAC_signal.bed',sep='\t',index=None)
    return None

def revise4():
    test = pd.read_csv('cluster_with_promoter.bed', sep='\t').fillna('0')
    group = test.groupby('start')

    output = pd.DataFrame(columns=['#chr','start','end','cluster','gene_id'])
    for sub_group in group:
        test_name = sub_group[0]
        df_test = sub_group[1]
        #print(df_test)
        for i in range(len(df_test)):
            if i >0: 
                df_test.iloc[0,4] = df_test.iloc[0,4] + ',' + df_test.iloc[i,4]
        output = pd.concat([output, df_test], axis=0)#output.append(df_test)

    output.drop_duplicates(subset=['#chr','start','end'],keep='first',inplace=True)
    output.to_csv('cluster_promoter_final.bed', sep='\t', index=None) 
    return None

def revise5():
    DNase = pybedtools.BedTool('ATAC_signal.bed')
    cluster = pybedtools.BedTool('cluster_promoter_final.bed')

    c1 = DNase.intersect(cluster, wa=True, wb=True)

    c1.saveas('D_C.bed')

    p = pd.read_csv('D_C.bed', sep='\t', names=['#chr_x','start_x','end_x','signal','#chr_y','start_y','end_y','cluster_y','gene_id_t'])
    p.drop_duplicates(subset=['#chr_y','start_y','end_y'],keep='first',inplace=True)
    c = pd.read_csv('cluster_promoter_final.bed', sep='\t')

    p = p.loc[:,['#chr_y','start_y','end_y','cluster_y','gene_id_t','signal']]

    c = c.merge(p, left_on='start', right_on='start_y', how='outer')

    c = c.loc[:,['#chr','start','end','cluster','gene_id','signal']]

    c = c.fillna(0)

    c.to_csv('cluster_annotated.bed', sep='\t', index=None)
    return c

def revise(col_fore, attach_name, promoter_name, PATH):
	test = pd.read_csv('cluster_'+col_fore+'.txt', sep='\t',names=['#chr','bins','cluster'], dtype={'#chr':str})
    #test = test[test['#chr'] == '1']
	output = pd.DataFrame(columns=['#chr','start','end','cluster','gene_id','signal'])
	group = test.groupby('#chr')

	for sub_group in group:
		df_test = sub_group[1]
		df_test.to_csv('cluster.txt',sep='\t',index=None,header=None)
		revise1(attach_name, PATH)
		revise2(PATH, promoter_name)
		revise3()
		revise4()
		x = revise5()
		output = pd.concat([output, x], axis=0)#output.append(x,ignore_index=None)
	output.to_csv('cluster_annotated_'+col_fore+'.txt',sep='\t',index=None)
        
	os.remove('ATAC_combine.bed')
	os.remove('ATAC_signal.bed')
	os.remove('cluster.txt')
	os.remove('cluster_combine.bed')
	os.remove('cluster_annotated.bed')
	os.remove('cluster_promoter_final.bed')
	os.remove('cluster_with_promoter.bed')
	os.remove('D_C.bed')
	os.remove('P_C.bed')
	#os.remove('cluster_'+col_fore+'.txt')
	return None

def run(argv):
    
	PATH_INPUT=argv.PATH + '/'
	col_fore = argv.label.split(',')[0]
	col_back  = argv.label.split(',')[1]
	promoter_path = argv.promoter_path
	other_path = argv.other_path


	print (" ")
	print("Run main")
	print ("Here is the Summary of your input.")
	print ("Input path of cluster files is: %s" % PATH_INPUT)
	print ("First label is: %s" % col_fore)
	print ("Second label is: %s" % col_back)
	print ('Your promoter file is: %s' % promoter_path)
	if other_path != 'no_123':
		print ('Your have another file is: %s' % other_path)
	print ("End of Summary.")
	print (" ")
	


#### Main 
	e=0
	f=0
	 
	if promoter_path == 'no_123':
		c = pd.DataFrame(data={'#chr':[],'start':[],'end':[],'gene_id':[]})
		c.to_csv(PATH_INPUT + 'promoter_illusion.bed', sep='\t', index=None)
		promoter_path = 'promoter_illusion.bed'
		e = 1
	if other_path == 'no_123':
		d = pd.DataFrame(data={'#chr':[],'start':[],'end':[],col_fore:[],col_back:[],'logFC':[]})
		d.to_csv(PATH_INPUT + 'other_illusion.bed', sep='\t', index=None)
		other_path = 'other_illusion.bed'
		f = 1
	revise(col_fore, other_path, promoter_path, PATH_INPUT)
	revise(col_back, other_path, promoter_path, PATH_INPUT)
	
	if e == 1:
		os.remove(PATH_INPUT + 'promoter_illusion.bed')
		
	if f == 1:
		
		os.remove(PATH_INPUT + 'other_illusion.bed')

	
	print(" ")
	return None

def main(argv):
	
	desc = "Combine cluster information with gene promoters and other genome functions."
	
	parser = OptionParser(description=desc)
	
	parser.add_option("-i", "--in", action="store", type="string",
			dest="PATH", help="Path to Input cluster, pomoter and another factor.", metavar="<file>")
	
	parser.add_option("-l", "--label", action="store", type="string",
			dest="label", help="Name of labels you have used to call hubs. Delimiter ',', for example 'label1,label2'.", metavar="<str>")
  
	parser.add_option("-p", "--promoter", action="store", type="string",
			dest="promoter_path", help="File name of gene promoters.", metavar="<file>")
	
	parser.add_option("-f", "--function_file", action="store", type="string", default = 'no_123',
			dest="other_path", help="File name of other factors.", metavar="<file>")


	(opt, args) = parser.parse_args(argv)
	if len(argv) < 3:
		parser.print_help()
		sys.exit(1)
	
	## parameters
	PATH_INPUT=opt.PATH + '/'
	
	col_fore = opt.label.split(',')[0]
	col_back  = opt.label.split(',')[1]
	
	promoter_path = opt.promoter_path
	other_path = opt.other_path


	print (" ")
	print("Run main")
	print ("Here is the Summary of your input.")
	print ("Input path of cluster files is: %s" % PATH_INPUT)
	print ("First label is: %s" % col_fore)
	print ("Second label is: %s" % col_back)
	print ('Your promoter file is: %s' % promoter_path)
	if other_path != 'no_123':
		print ('Your have another file is: %s' % other_path)
	print ("End of Summary.")
	print (" ")
	


#### Main
	e=0
	f=0
	 
	if promoter_path == 'no_123':
		c = pd.DataFrame(data={'#chr':[],'start':[],'end':[],'gene_id':[]})
		c.to_csv(PATH_INPUT + 'promoter_illusion.bed', sep='\t', index=None)
		promoter_path = 'promoter_illusion.bed'
		e = 1
	if other_path == 'no_123':
		d = pd.DataFrame(data={'#chr':[],'start':[],'end':[],col_fore:[],col_back:[],'logFC':[]})
		d.to_csv(PATH_INPUT + 'other_illusion.bed', sep='\t', index=None)
		other_path = 'other_illusion.bed'
		f = 1
	revise(col_fore, other_path, promoter_path, PATH_INPUT)
	revise(col_back, other_path, promoter_path, PATH_INPUT)
	
	if e == 1:
		os.remove(PATH_INPUT + 'promoter_illusion.bed')
		
	if f == 1:
		
		os.remove(PATH_INPUT + 'other_illusion.bed')

	
	print(" ")
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)