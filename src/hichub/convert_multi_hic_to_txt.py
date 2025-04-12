import pandas as pd
import os
import struct
from optparse import OptionParser
import sys
#import straw
import hicstraw


def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return buf.encode("utf-8", errors="ignore")
            return buf.decode("utf-8")
        elif b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        else:
            buf += b
    return None

def read_header(req):
    """
    Takes in a .hic file and returns a dictionary containing information about
    the chromosome. Keys are chromosome index numbers (0 through # of chroms
    contained in file) and values are [chr idx (int), chr name (str), chrom
    length (str)]. 
    """
    chrs = {}
    resolutions = []
    magic_string = struct.unpack(b'<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        print('This does not appear to be a HiC file.')
        exit(0)
        #force_exit(error_string, req)
    global version
    version = struct.unpack(b'<i', req.read(4))[0]
    
    masterindex = struct.unpack(b'<q', req.read(8))[0]
    genome = b""
    c = req.read(1)
    while (c != b'\0'):
        genome += c
        c = req.read(1)
    genome = genome.decode('ascii')
    # metadata extraction
    metadata = {}
    nattributes = struct.unpack(b'<i', req.read(4))[0]
    for x in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
        metadata[key] = value
    nChrs = struct.unpack(b'<i', req.read(4))[0]
    for i in range(0, nChrs):
        name = readcstr(req)
        length = struct.unpack(b'<i', req.read(4))[0]
        if name and length:
            chrs[i] = [i, name, length]
    nBpRes = struct.unpack(b'<i', req.read(4))[0]
    # find bp delimited resolutions supported by the hic file
    for x in range(0, nBpRes):
        res = struct.unpack(b'<i', req.read(4))[0]
        resolutions.append(res)
    return chrs, resolutions, metadata
    
def HiC_Matrix_to_Txt(_NORM, _hic, _cond, _resolution, _chr, _cut_off):
    result = hicstraw.straw('observed', _NORM, _hic, _chr, _chr, 'BP', _resolution)
    x=[]
    y=[]
    c=[]
    for i in range(len(result)):
        x.append(result[i].binX)
        y.append(result[i].binY)
        c.append(result[i].counts)
    df_tem = pd.DataFrame(data={'bin1':x,'bin2':y, _cond:c})  
    df_tem = df_tem[df_tem.loc[:, _cond] > _cut_off]
    #df_tem.loc[:, _cond] = round(10**6*df_tem.loc[:, _cond] /  df_tem[(df_tem.bin1==df_tem.bin2)].loc[:, _cond].sum(),2)
    return df_tem.fillna(0)
    
def Multi_Input_Matrix_to_Txt(_Norm, _hics, _conds, _resolution):
    
    Out_Name = 'Summary_'+'_'.join(_conds)+'_Dense_Matrix.txt'
    df_hic = pd.DataFrame(columns=['#chr','bin1','bin2']+_conds).to_csv(Out_Name, sep='\t',
                                                                      header=True, index=None)
    req = open(_hics[0], mode='rb')
    chrs, resolutions, metadata = read_header(req)
    for idx in chrs:
        chr_idx = str(chrs[idx][1])#.replace('chr','')
        if(chr_idx not in {'ALL','All', 'chrM','M'}):
            df_hic = pd.DataFrame(columns=['bin1','bin2'])
            for cond, hic in zip(_conds, _hics):
                df_hic_tem = HiC_Matrix_to_Txt(_Norm, hic, cond, _resolution, chr_idx, 0)
                if len(df_hic_tem) != 0:
                    df_hic = df_hic.merge(df_hic_tem, on=['bin1','bin2'], how='outer')
                else:
                    df_hic_tem = pd.DataFrame(columns=['bin1', 'bin2', cond])
                    df_hic = df_hic.merge(df_hic_tem, on=['bin1','bin2'], how='outer')
                
            df_hic.insert(loc=0, column='#chr', value=chr_idx)
            df_hic = df_hic.fillna(0)
            df_hic = df_hic.loc[:,['#chr','bin1','bin2']+_conds]
            df_hic.to_csv(Out_Name, sep='\t', mode='a', header=False, index=None)
    return None

### FUNCTION
####################################################################################
### FUNCTION
### FUNCTIONS
def run(opt):
	
	PATH_INPUT = opt.input_path
	os.chdir(PATH_INPUT)
	File_Name_Sets = opt.file_name.split(',')
	File_Label_Sets = opt.file_label.split(',')
	resolution = opt.res
	
	Norm_Meths = opt.norm_hic  ## "KR" or "NONE"
	
	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input path of HiC files in .hic format: %s" % opt.input_path)
	print ("Name of input .hic files: %s" % opt.file_name)
	print ("Label of output .txt files: %s" % opt.file_label)
	#print ("Norm of Input File: %s" % opt.norm_hic)
	print ("Resolution is: %i" % opt.res)
	print ("End of Summary.")
	print (" ")

#### Main 
	Multi_Input_Matrix_to_Txt(Norm_Meths, File_Name_Sets, File_Label_Sets, resolution)
	return None

def main(argv):
	
	desc = "Convert multi .hic to .txt format."
	
	parser = OptionParser(description=desc)
	
	parser.add_option("-i", "--input_path", action="store", type="string",
			dest="input_path", help="Path to input HiC files.", metavar="<path>")
	
	parser.add_option("-f", "--file_name", action="store", type="string",
			dest="file_name", help="Name of input .hic files. Delimiter ',', for example 'file1,file2'.", metavar="<files>")
	
	parser.add_option("-l", "--file_label", action="store", type="string",
			dest="file_label", help="Label of output .txt files. Delimiter ',', for example 'label1,label2'.", metavar="<str>")
	
	parser.add_option("-r", "--resolution", action="store", type="int",
		dest="res", help="Resolution of output contact matrixs.", metavar="<int>")
	
	
	parser.add_option("-n", "--norm", action="store", type="string", default = 'NONE',
			dest="norm_hic", help="Normalization method of reading .hic files, default=NONE.", metavar="<str>")

	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input path of HiC files in .hic format: %s" % opt.input_path)
	print ("Name of input .hic files: %s" % opt.file_name)
	print ("Label of output .txt files: %s" % opt.file_label)
	#print ("Norm of Input File: %s" % opt.norm_hic)
	print ("Resolution is: %i" % opt.res)
	print ("End of Summary.")
	print (" ")
	
	## parameters
	PATH_INPUT = opt.input_path
	os.chdir(PATH_INPUT)
	File_Name_Sets = opt.file_name.split(',')
	File_Label_Sets = opt.file_label.split(',')
	Norm_Meths = opt.norm_hic  ## "KR" or "NONE"
	resolution = opt.res

#### Main 
	Multi_Input_Matrix_to_Txt(Norm_Meths, File_Name_Sets, File_Label_Sets, resolution)
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
