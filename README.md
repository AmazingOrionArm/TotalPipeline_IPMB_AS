# DNA Methylation Toolkit
DNA Methylation Toolkit provide reseachers to find out the relationship of methylation in different samples.
There are many tools contain is this Toolkit, users can use it discretionary.

# Tools in Toolkit
* Heatmap & PCA：
* Differentially Methylated Region：
* Fold Enrichment：
* Differentially Methylated Gene：
* Chromosome View：
* Metaplot：

# Before using Toolkit
Before using Toolkit, system requirement and some kits that need to be installed are as follows
* System Requirement
  * CPU：No special restrictions, but CPU has 16 cores is more efficient
  * MEM：12GB or higer (in arabidopsis sample) / 256GB or higher (in human sample)

* Essential Kits
  * Python **2.7** or above
  * Python modules：pandas, numpy, matplotlib, math, scipy, argparse, glob, pyBigWig, collections, gzip, re, PyQt5
  * R **3.5.1** or above
  * bedtools
  * deepTools
  
## Download and install python
for windows
>     https://www.python.org/downloads/%20windows/
  * and choose version which environment you are

for Mac OS X
>     https://www.python.org/downloads/mac-osx/
  * and choose version which environment you are

for Linux (Ubuntu, also in Mac OS X), using command line
>     sudo apt-get install python2.7
or
>     sudo apt-get install python3.6

## Download and install R
In this link, choose one mirror that is near your location, and choose version which environment you are
>     https://cran.r-project.org/mirrors.html

or in Mac OS X and Linux, using this command
>     sudo apt-get install r-base

## Download and install pip
for windows
>     python get-pip.py

for Mac OS X
>     sudo apt-get install python-pip

for Linux (ubuntu)
>     sudo apt-get install python-pip

## Download and install essential kits
Use the following commands one by one
>     pip install pandas
>     pip install numpy
>     pip install matplotlib
>     pip install math
>     pip install scipy
>     pip install argparse
>     pip install glob
>     pip install pyBigWig
>     pip install collections
>     pip install gzip
>     pip install re
>     pip install PyQt5

# How to use
User can choose nongui or gui mode to use this Toolkit
* GUI mode
  * See GUI description.docx for details.  

* Non-gui mode
using commandline
  * program path：Path of 'nongui.py' (This Toolkits main program(non-gui mode))
  * input_gene_name：Path of input_gene (Your input gene file)
  * Optional Parameters：Analysis parameters. If you not use it, it will setting by default value
>     python [program path] [input_gene_name] [Optional Parameters]

Parameters of this program
>     positional arguments:
>       input_gene_name       path of input gene

>     optional arguments:
>       -h, --help            show this help message and exit
>       -d DEPTH              min size of #C+#T
>       -r REGION             size of region
>       -q QUAIFIED           min number of region member
>       -hc HEATMAP_CUTOFF    Heatmap_cutoff
>       -dcgc DMR_CG_CUTOFF   DMR_CG_cutoff
>       -dchhc DMR_CHH_CUTOFF DMR_CHH_cutoff
>       -dchgc DMR_CHG_CUTOFF DMR_CHG_cutoff
>       -b BIN_SIZE           resolution of chrView and Metaplot
>       -p PROMOTER_SIZE      promoter_size

For example, if a user want to analysis arabidopsis sample, and he want to set region size and dmr's cut off value, he need input this command
>     python nongui.py arabidopsis.gtf -r 800 -dchhc 0.1 -dchgc 0.02
