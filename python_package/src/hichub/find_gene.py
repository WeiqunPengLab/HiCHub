import pandas as pd
import numpy as np
import igraph as ig
#import pybedtools
from scipy import stats
from optparse import OptionParser
import sys, os, multiprocessing
import random
pd.options.mode.chained_assignment = None  # default='warn'

def revise_hub(hubs):
    df_test = hubs
    #df_test = df_test[df_test['-log10(pvalue)']>20]
    pyramid = pd.DataFrame()
    pyramid = df_test[df_test['reg1'] == df_test['reg2']]

    pyramid.loc[:,'0'] = pyramid.loc[:,'reg1'].str.split(':', expand = True)[0]
    pyramid.loc[:,'1'] = pyramid.loc[:,'reg1'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    pyramid.loc[:,'2'] = pyramid.loc[:,'reg1'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    stripe_1 = pd.DataFrame()
    stripe_1 = df_test[df_test['reg1'] != df_test['reg2']]
    stripe_1.loc[:,'0'] = stripe_1.loc[:,'reg1'].str.split(':', expand = True)[0]
    stripe_1.loc[:,'1'] = stripe_1.loc[:,'reg1'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    stripe_1.loc[:,'2'] = stripe_1.loc[:,'reg1'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)

    stripe_2 = pd.DataFrame()
    stripe_2 = df_test[df_test['reg1'] != df_test['reg2']]
    stripe_2.loc[:,'0'] = stripe_2.loc[:,'reg2'].str.split(':', expand = True)[0]
    stripe_2.loc[:,'1'] = stripe_2.loc[:,'reg2'].str.split(':', expand=True)[1].str.split('-', expand=True)[0].astype(int)
    stripe_2.loc[:,'2'] = stripe_2.loc[:,'reg2'].str.split(':', expand=True)[1].str.split('-', expand=True)[1].astype(int)


    stripe = pd.DataFrame()
    stripe = stripe.append(stripe_1)
    stripe = stripe.append(stripe_2)
    stripe.drop_duplicates(subset=['0','1','2'],keep='first', inplace=True)


    total = pd.DataFrame()
    total = total.append(stripe)
    total = total.append(pyramid)
    total.drop_duplicates(subset=['0','1','2'],keep='first', inplace=True)
    total['label'] = 'H1ESC'

    total = total.loc[:,['0','1','2','reg1','reg2','-log10(pvalue)','label']]
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
        sys.exit(1)
    
    k=0
    for i in range(len(file_name)):
        hub = pd.read_csv(input_path+file_name[i], sep='\t')
        hub = revise_hub(hub)
        hub = hub[hub['0'] == chrom]
        hub = hub[hub['1'] <= start]
        hub = hub[hub['2'] >= start]
        hub = hub.loc[:,['reg1','reg2','-log10(pvalue)']]
        if len(hub) > 0:
            print('The gene ' + gene_name + ' is in '+file_label[i]+' hubs:')
            k = k+1
            x=[]
            for i in range(len(hub)):
                x.append(i+1)
            print(hub.to_string(index=False))
    if k == 0:
        print('The '+gene_name+' does not exist in any hubs.')
    return None
    

####################################################################################
### FUNCTION
### FUNCTIONS
def run(opt):
	## parameters
	PATH_INPUT = opt.input_path
	os.chdir(PATH_INPUT)
	File_Name_Sets = opt.file_name.split(',')
	File_Label_Sets = opt.file_label.split(',')
	gene_name = opt.gene_name  ## "KR" or "NONE"

	
	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in .hic format: %s" % opt.input_path)
	print ("Name of Input File: %s" % opt.file_name)
	print ("Label of Input File: %s" % opt.file_label)
	print ("Name of gene: %s" % opt.gene_name)
	print ("End of Summary.")
	print (" ")
	
#### Main 
	find_hub(gene_name, PATH_INPUT, File_Name_Sets, File_Label_Sets)

	print(" ")
	return None

def main(argv):
	desc="Please enter the path of Hichub cluster output file, converted file form hic and gene name."
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Path to Input HiC file in txt format", metavar="<file>")
	parser.add_option("-n", "--gene", action="store", type="string",
			dest="gene_name", help="Name of gene", metavar="<str>")
	parser.add_option("-f", "--file_name", action="store", type="string",
			dest="file_name", help="Name of hub's File. Delimiter ',', for example 'file1,file2' ", metavar="<str>")
	parser.add_option("-l", "--file_label", action="store", type="string",
			dest="file_label", help="Label of hub. Delimiter ',', for example 'label1,label2'", metavar="<str>")
    
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	## parameters
	PATH_INPUT = opt.input_path
	os.chdir(PATH_INPUT)
	File_Name_Sets = opt.file_name.split(',')
	File_Label_Sets = opt.file_label.split(',')
	gene_name = opt.gene_name  ## "KR" or "NONE"



	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in .hic format: %s" % opt.input_path)
	print ("Name of Input File: %s" % opt.file_name)
	print ("Label of Input File: %s" % opt.file_label)
	print ("Name of gene: %s" % opt.gene_name)
	print ("End of Summary.")
	print (" ")
	


#### Main 
	find_hub(gene_name, PATH_INPUT, File_Name_Sets, File_Label_Sets)

	
	print(" ")
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
