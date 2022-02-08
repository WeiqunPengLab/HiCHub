# NOTE: For Paper Review, please follow the instruction below:
Latest updated on 01/21/2022,

# Comprehensive Network Analysis for HiC

<br>
<img src="image/Hub_Myb.png" width="800">
<br>

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Example of Running (Demo)](#Example_Running)
- [License](#license)

## Overview
This module is a Python package containing tool for network analysis of HiC data.
It starts from HiC Interaction pairs, then generating network and clustering. Finally ranking all clusters by their interaction change.

## System Requirements
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
```

## Installation Guide

```
Quick install: python -m pip install git+https://github.com/WeiqunPengLab/hic_hub_test
```

Alternative method of installation:
Recommend to use bioconda for installing.
Create the environment from the environment_hichub.yml (Can be found in this repository) file:
```
conda env create -f environment_hichub.yml
python3 -m pip install hichub --user
python3 -m pip install numpy pandas pybedtools python-igraph scipy
```

```
https://bioconda.github.io/user/install.html
```

## Example of Running (Demo)
After installation, type hichub in your command line will return the following UI:
```
welcome
The python env is: 3.7.10 | packaged by conda-forge | (default, Feb 19 2021, 16:07:37)
[GCC 9.3.0]
usage: hichub [-h] [-v] {test,diff,convert} ...
hichub -- A toolset to detect and analyze differential Hubs
positional arguments:
  {test,diff,convert}  sub-command help
    test               Output for test parser
    convert            Convert multi .hic to txt format --> Format should be: #chr bin1 bin2 Count
    diff               Parser for diff hubs
    find               Find gene in which hubs
    plot               Draw community plot to find gene in hubs
optional arguments:
  -h, --help           show this help message and exit
  -v, --version        show program's version number and exit
```

diff:
Find hichubs for different cell conditions of the same specious.
```
usage: hichub diff [-h] -i <file> -f <str> -b <str> -r <int> [-d <float>] [-c <float>] [-p <float>] [-t <int>]

-i: Path to Input HiC file in txt format.
-f: Name of condition as foreground.
-b: Name of condition as background.
-r: Resolution if Hic .txt
-d: Foldchange to evaluate difference, default 1.0.
-c: Cut-off value to remove low abundance interaction pairs, default 10.0.
-p: Pvalue cutoff for output (diff hub), default 1e-5.
-t: Number of threads to run, default 1.

hichub diff --help 
Input Format: HiC Interaction in txt format. Example of test data can be found in ~/test_data
And the output can be found at working directory: 
(Demo output is: 6_1.0_DKO_na_WT_na_specific_regions.bed or 4_1.0_WT_na_DKO_na_specific_regions.bed)
```

convert:
Convert (.hic) files to required input format
```
usage: hichub convert [-h] -i <file> [-n <str>] -f <str> -l <str> -r <int>

Exmaple of output: 
reg1    reg2    -log10(pvalue)
chr10:19570000-19850000 chr10:19570000-19850000 45.91186483627381
chr10:20860000-21060000 chr10:20860000-21060000 41.129022601906215
chr10:117030000-117140000       chr10:116870000-117010000       14.165
chr10:95130000-95290000 chr10:95130000-95290000 9.80623027538454
chr10:18970000-19160000 chr10:18970000-19160000 9.288099829570816
```

find:
Find gene in which hub.
```
usage: hichub find [-h] -i <input file path> -f <str> -l <str> -n <gene name>
-i: Path to store all required input files.
-f: File names for different hubs, the output files for 'diff' call hubs. <file_name_1,file_name_2,...>
-l: The label for your input files. <condition_1, condition_2,...>
-n：The gene's name that you want to find in hubs. Examples: for human, LSR, SOX2, SERPINB8
```

plot:
Plot community igraph-plot for the community that contain the gene you want to find in your hubs.
```
usage: hichub plot [-h] -i <input path for your .hic file> -c <input cluster files form 'diff'> -n <gene name> -f <str> -b <str>
-i: Path to Input HiC file in txt format.
-f: Name of condition as foreground.
-b: Name of condition as background.
-n: Gene name that you want to plot.
-c: Path to Input cluster file calcuated from 'diff'.
```

## Built With

## Contributing

Please read (https:xx) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

## Authors

* *Xiang Li  

## Contributors
* Shuang Yuan
## License

#This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments


