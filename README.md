Latest updated on 29/03/2022,

# Comprehensive Network Analysis for HiC

<br>
<img src="image/Hub_Myb.png" width="800">
<br>

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Example of Running (Demo)](#Example_Running)
- [License](#license)

## Overview
This module is a Python package containing tool for network analysis of HiC data.
It starts from HiC Interaction pairs, then generating network and clustering. Finally ranking all clusters by their interaction change.

### Hardware Requirements

This package requires only a standard computer with enough RAM to support the in-memory operations.

### Software Requirements

HicHub mainly depends on the Python scientific stack.

```
python >=3
pandas
numpy
pybedtools
python-igraph
scipy
hic-straw
statsmodels
```

## Installation Guide

```
Quick install: python -m pip install git+https://github.com/WeiqunPengLab/HicHub_installation
```

## Example of Running (Demo)

This tool contains 5 functions, which are 'convert', 'diff', 'asso', 'find', 'plot'. The examples of using are as follows:

1) 'convert'
  Please prepare two '.hic' files for comparing. Put these two files in the same directory. Then run following code:
  
```
hichub convert -i [the directory of your files] -f [XX1.hic,XX2.hic](your files name) -l [XX1,XX2](name your file) -r [resolution] -n [normalization method]
```

  ** -n normaliztion method is optional, the default is 'NONE', you can try 'KR'. The details can refert to <hic-straw>.
 
2) 'diff'
  Please put your converted file in setp (1), such as 'chr1.txt' (test data) in a directory and run in the same directory as:

```
hichub diff -i [the file name you have in step (1)] -f [one name of your label] -b [the other name of your label] -r [resolution]
```
  ** Here are some optional parameters:
  -p p-value you want to cut, default = 0.00001
  -d Fold change of difference of two situation to do cluster analysis, default = 1.0
  -c Cut-off value, to remove the value of the sum of your two situation smaller than, default  = 10
  
  The example for test data:

```
hichub diff -i chr1.txt -f H1ESC -b HFFc6 -r 10000 -c 10 -d 1 -p 0.000001
```
  
3) 'associate'
  After run step (2), you will get two files name as "cluster_XX1.txt", "cluster_XX2.txt", them under this directory, give the promoter file and a function file (such as DNase, CTCF...), make them in the following format
  Promoter file: #chr (Tab) start (Tab) end (Tab) gene_id
  Other function file: #chr (Tab) start (Tab) end (Tab) XX1 (Tab) XX2 (Tab) logFC--[log2(XX1/XX2)]
  There are two example files in the test folder under this github, see the format of them.
  
  Then run:
  
```
hichub asso -i [the directory you store these files)] -f [one name of your label] -b [the other name of your label] -p [promoter files name] -o [function files name]
```
  
  The example for test data:

```
hichub asso -i [] -f H1ESC -b HFFc6 -p promoter.bed -o DNase.bed
```  
  
4) 'find'
  After run setp (2), if you have a gene promoter file, you can search gene in out hubs:

```
hichub find -i [the directory you store resultes of step 2 and promoter file)] -f [your hub results name XXX1.bed,XXX2.bed] -l [label your results] -n [name of your gene]
```
  
5) 'plot'
  If you run step (3) before, you can use two cluster final files to plot a igraph-plot for your hub, you need the chrmome number 'chr?' and the first number of your hubs, such as output 'chr3:104700000-105380000', you need 'chr3' and '104700000'
  
  Put all your files in the same directionary and run:

```
hichub plot -i [your original .txt file] -f [XX1] -b [XXX2] -d [FC you us before] -c [cut_off value you use before] -o [your culseter file] -p [chrom, such as 'chr3'] -n [node, such as 104700000]
```
optional -g [if you have gene file and enter a gene name in step (3)]  


## Versioning

```
1.0.0
```

## Authors

* Xiang Li, Shuang Yuan, Shaoqi Zhu  

## License

#This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
