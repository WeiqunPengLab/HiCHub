Latest updated on 01/10/2022,

# Comprehensive Network Analysis for HiC

<br>
<img src="image/Hub_Myb.PNG" width="800">
<br>

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Example of Running](#Example_Running)
- [License](#license)

## Overview
This module is a Python package containing tool for network analysis of differential interactions of Hic data.  

### Hardware Requirements

This package requires only a standard computer with enough RAM to support the in-memory operations.

### Software Requirements

We strong recommand you to install Anaconda before install HicHub. 
HicHub mainly depends on the Python scientific stack.  

```
python >=3
pandas = 1.4.3
numpy = 1.23.0
bedtool >= 2.70.1
pybedtools = 0.9.0
python-igraph = 0.9.11
scipy = 1.8.1
hic-straw = 1.3.1
statsmodels = 0.13.2
pycairo >= 1.11.0
```

## Installation Guide
Before install HicHub, make sure you install following packages:
```
pip install hic-straw
```
```
sudo apt-get install bedtools
```
```
pip install pybedtools
```
```
pip install pycairo
```
```
pip install scipy
```

Quick installation. Please type the following command in your Linux shell to install HicHub.

```
python -m pip install git+https://github.com/WeiqunPengLab/HiCHub
```

If installing successfully, type

```
hichub
```

in your Linux shell, then you will see the following interface:

```
welcome
The python env is: [your python version] (main, [time you type this command])
[GCC 7.5.0]
usage: hichub [-h] {convert, diff, asso, plot} ...
hichub -- A toolset to detect and analyze differential Hubs.
positional arguments:
  {convert, diff, asso, plot}  sub-command help
    convert            Convert multi .hic to .txt format.
    diff               Parser for call hubs
    asso               Parser for associate clusters with promoter
    plot               Parser for plot clusters in igraph package.
optional arguments:
  -h, --help           show this help message and exit
```
  
## Example of Runing

In order to call HicHubs, you first need to prepare two (.hic) files and put them in the same directory.  
  
In this Github, under the directory '~/test' there are two (.hic) files named 'H1ESC.hic' and 'HFFc6.hic'.  
  
Please download them for test.  
  
ATTENTION: You need to run all these processes in the same directory, once you finished one step, please don't change the name of output files!  
  
### 'convert'
First, you need to use 'convert' command to convert them to (.txt) format that we will use,

```
hichub convert -i [run path] -f [file names, seperate with ','] -l [label of output, seperate with ','] -r [resolution of bin]
```
 -i : Your input path, the directory that you strore your two (.hic) files, the directory that you run the HicHub program.  
 -f : Your input file names, seperate them with comma. For example '-f H1ESC.hic,HFFc6.hic'.  
 -l : Your output files' labels, name your two input files with the same order, seperate them with comma. For example '-l H1ESC,HFFc6'.  
 -r : The length for one bin on genome , the unit is 'bp'. For example '-r 10000'.  
 -n : (optional, default='NONE'), the normalizaition method for converting .hic files, consistent with the method provided by hic-straw. 'NONE' for no normalization and 'KR' for KR normalization.

For example: (The test data may takes you around 2 mins to finish.)  

```
hichub convert -i /mnt/d/test -f H1ESC.hic,HFFc6.hic -l H1ESC,HFFc6 -r 10000
```

The output is a (.txt) format files which contains the contact matrics of two (.hic) files in the format:  

```
#chr    bin1    bin2    label1    label2
```
(The empty between two labels is tab.)  
  
For example (the output of test data):  
```
#chr    bin1    bin2    H1ESC    HFFc6  
```

Where, '#chr', 'bin1', 'bin2' represent the chromosome, the location of the left anchor and the location of the right anchor respectivly. 
  
<img src="image/3.png" width="800">
  
### 'diff'  
When you converted two (.hic) files into the (.txt) format by command 'convert', you can use command 'diff' to call hubs.

```
hichub diff -i [yout (.txt) file's name] -l [label you have used before] -r [resolution of bin] -c [optional, cut-off threshold] -d [optional, folde change threshold] -p [optional, p-value threshold] -n [optional, normalization method]
```
 -i : Your converted (.txt) file's name from the 'convert' step. For example: '-i Summary_H1ESC_HFFc6_Dense_Matrix.txt'  
 -l : The labe you have used in the 'convert' step. For example: '-l H1ESC,HFFc6'  
 -r : The resolution you have used in the 'convert' step. For example: '-r 10000'  
 -n : Optional default = 'LOESS', please choose normalization method from 'LOESS' or 'Dis_LOESS'(doing LOESS normalizaiton at each distance). The details could be   seen in the paper.  
 -c : Optional default = 10, remove the sum of two values of contact matric samller than a threshold. For example: '-c 10'.  
 -d : Optional default = 0, the difference to determine whether we keep an edge in cluster analysis. The details could be seen in the paper.  
 -p : Optional default = 0.00001, The threshold to pick hubs smaller than a certain p-value.  

For example: (The test data may takes you around 4 mins to finish.)  
```
hichub diff -i Summary_H1ESC_HFFc6_Dense_Matrix.txt -l H1ESC,HFFc6 -r 10000
```
```
hichub diff -i Summary_H1ESC_HFFc6_Dense_Matrix.txt -l H1ESC,HFFc6 -n Dis_LOESS -r 10000 -c 10 -d 0 -p 0.00001
```
There are four output files.    
(1) --- 'H1ESC_specific_hubs_comparing_with_HFFc6.bed'  
(2) --- 'HFFc6_specific_hubs_comparing_with_H1ESC.bed'  
(3) --- 'cluster_H1ESC.txt'  
(4) --- 'cluster_HFFc6.txt'

(1) and (2) contain the cell-type-specific hubs we found. The format of output file (hubs) is:  
```
#chr    left_hub_anchor_start    left_hub_anchor_end    right_hub_anchor_start    right_hub_anchor_end    -log10(pvalue)    -log10(FDR)  
```
(The empty between two labels is tab.)  
  
(3) and (4) record the cluster information we used to call hubs, they are use for drawing the network plot in the following functions.  
  
<img src="image/4.png" width="800">
  
### 'asso'
This function associate genome annotation (gene, open chromatin) with network clusters. 
  
In the '~/test' folder, there are also two test files named 'promoter.bed' and 'DNase.bed' containing the information of gene promoter and DNase for the test data.  
  
```
hichub  asso -i [run path] -l [label you have used before] -p [the files contain gene promoter] -f [Optional, file name for DNase, CTCF, ...]
```
 -i : Your input path,the directory that you run the HicHub program. For example: '-i /mnt/d/test'  
 -l : The labe you have used in the 'convert' and 'diff' step. For example: '-l H1ESC,HFFc6'  
 -p : The name of the file that contains gene promoter's information. It contains the coordinate and name of the genes you want to analize and plot in the future. The input format should be :
```
#chr    start    end    gene_name
```     
(The empty between two labels is tab.)  
  
 -f : Optional. The name of the file that contains information about another factor, such as DNase, CTCF ....
      The input format shows as following, where signal represents if your factor's count is increasing or decreasing between two situation, if in situation label1 larger than label2, please mark them as 'up', if in situation label2 larger than label1, please mark them as 'down', if not change, please mark them as number 0. The label1 and label2 represent the counts of your factor in the given regions.   
```
#chr    start    end    label1    label2     signal  
```  
(The empty between two labels is tab.)  
 
For example: (The test data may takes you around 1 mins to finish.)  
```
hichub asso -i /mnt/d/test -l H1ESC,HFFc6 -p promoter.bed
```
```
hichub asso -i /mnt/d/test -l H1ESC,HFFc6 -p promoter.bed -f DNase.bed
```
The ouput is the files with name "cluster_annotated_H1ESC.txt" and "cluster_annotated_H1ESC.txt". They associate the gene with the cluster nodes.
```
The first three columns of the output files are the coordinates of each node. The forth column ("cluster") marks which cluster this node belongs to. The fifth column ("gene_id") records if there are any gene include in this node. If there are no gene in this node, it will record the "0". If there are multi-gene in this node, their names will be seperated by ",". The sixth column "signal" depands on which kind of input file you choose for "-f" option. It will give you "up" (for increasing the expresstion), "down" (for decreasing the expresstion) or "0" (for does not change).
  
<img src="image/5.png" width="800">
  
### 'plot'
This function plots the networks associated with one or more particular gene(s) along with the annotation information.

```
hichub plot -i [yout (.txt) file's name] -l [label you have used before] -p [the files contain gene promoter] -n [gene names, seperate with ','] -c [optional, cut-off threshold you have used in 'diff'] -d [optional, folde change threshold you have used in 'diff'] 
```
 -i : Your converted (.txt) file's name from the 'convert' step. For example: '-i Summary_H1ESC_HFFc6_Dense_Matrix.txt'  
 -l : The labe you have used in the 'convert', 'diff', 'asso' steps. For example: '-l H1ESC,HFFc6'  
 -p : The file name that contain gene promoter's information.    
      The input format should be : #chr----start----end----gene_name   
 -n : The gene names you want to plot. For example: '-n CPTC,DFFB'  
 
 -c : Optional default = 10. The -c value you have used in 'diff'.       
 -d : Optional default = 1.0. The -d value you have used in 'diff'.     
 
The result will contain the network plot for gene you entered.   
If this gene is not in any hub, you will received the hint that : 'The 'gene_name' does not exist in any hubs.'  
For example: (The test data may takes you around 3 mins to finish.)  
```
hichub plot -i Summary_H1ESC_HFFc6_Dense_Matrix.txt -l H1ESC,HFFc6 -p promoter.bed -n CPTP
```
It will plot the network around gene CPTP.  

<img src="image/6.png" width="800">
## Versioning

```
2.0.0
```

## Authors

* Xiang Li, Shuang Yuan, Shaoqi Zhu  

## License

#This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
