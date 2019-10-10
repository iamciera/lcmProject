# Preprocessing RNAseqAnalysis_LCM

This my lab notebook on my RNAseq analysis for my Laser Capture Micro-dissection (LCM) project.

License: [Attribution-NonCommercial-ShareAlike 4.0 International](http://creativecommons.org/licenses/by-nc-sa/4.0/)
Author: Ciera Martinez

## Overview of preprocessing of Illumina Reads

Each of the three lanes were processed seperately.  The preprocessing of raw reads was done on the [iplant atmoshere](http://www.iplantcollaborative.org/ci/atmosphere) image Maloof08 v2.  The full analysis was then stored on [IRODS](http://irods.org/).  

Steps performed on each lane

1.  Checkpoint 1: `fastqc` to understand overall quantity and quality of raw reads.
2.  `trimFastQuality.py` to trim sequences and reject based on low quality. 
3.  `read_N_remover.py` to remove sequences containg "N" nucleotides. 
4.  `adapterEffectRemover.py` remove adapter contamination.
5.  Checkpoint 2: `fastqc` to understand what the preprocessing did to reads. 
6.  `barcode_split_trim.pl` 
    a. splits reads by barcode
    b. sorts each read into seperate fasta file based on barcodes
    c. gives read count per barcode

# Detailed protocol for Illumina data preprocessing for in-line barcoded data 

## Packages

[http://comailab.genomecenter.ucdavis.edu/index.php/Barcoded_data_preparation_tools](http://comailab.genomecenter.ucdavis.edu/index.php/Barcoded_data_preparation_tools)

Information about the settings and usage can be found in the text file of the scripts themselves.  

### fastQC 

Already installed on the Maloof instance.  Documention can be found at [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt).

### fastx



### auto_barcode

On Github: `git clone https://github.com/mfcovington/auto_barcode.git`

#### Installing

You will also require [fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html), from the [Hannon Lab](http://hannonlab.cshl.edu/). For the iPlant Atmosphere environment download the 64 bit linux version. Pre-compiled binaries available for: Linux (64 bit).  Followed directions on page and put the folder in my scripts folder. 

[pkg-config guide](http://people.freedesktop.org/~dbn/pkg-config-guide.html)
The primary use of pkg-config is to provide the necessary details for compiling and linking a program to a library.

Download pre-compiled binaries, put them in /usr/local/bin
	
	$ mkdir fastx_bin
	$ cd fastx_bin
	$ wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
	$ tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
	$ sudo cp ./bin/* /usr/local/bin
	

## Preprocessing Protocol

The files come in .gz format so first you need to unzip them in your working folder.  `gunzip` replaces the original zipped file with the unzipped to save space.

	$ gunzip *.gz #10 min

The files then need to be concatenated into a single file  .fq or .fastq.

	$ cat *.fq > MyLane.fq #10 min

To get sequence quality we ran FastQC and to compare total vs. filtered reads. FastQC puts together information on reads in a tidy .html file with graphs.  Later we will did the same after quality filtering to get an understanding of what sort of effects quality filtering had on our reads. 

	fastqc MyLane.fq  #~20 min

Trim sequences and reject based on low quality scores.

	$ python trimFastqQuality.py 20 35 MyLane.fq QualFiltered.fq # Time takes around ~2 hours per lane

We then ran FastQC on `QualFiltered.fq` and compared the total reads to the  unfiltered read file.

We removed sequences containing “N” nucleotides

	$ python read_N_remover.py QualFiltered.fq  Nremoved.fq #15 minutes

and remove adapter contamination sequences.

	$ python adapterEffectRemover.py 41 Nremoved.fq AdaptersRemoved.fq b #10 min

FastQC on `AdaptersRemoved.fq` file.

	fastqc AdaptersRemoved.fq #~20 min

### Barcode Splitting

Fastx:

	$ cat AdaptersRemoved.fq | perl fastx_barcode_splitter.pl --bcfile BCfile.txt --bol --exact --prefix Split_ --suffix ".fq" # time = ~ > 5 hours

### Trim Barcodes using auto_barcode

You first need to switch the  columns of your BCfile.txt, because auto_barcode needs them a different way. Ie.

ATAGG	Barcode1	
GCTAT	Barcode2	

	awk '{ print $2 "\t" $1}' BCfiletest.txt > BCfile2.txt #switches the columns

To get the scripts

    git clone https://github.com/mfcovington/auto_barcode.git

Download Perl dependency

	cpanm --sudo Text/Table.pm

Run code 

    barcode_split_trim.pl [options] --barcode barcode.file IN.FASTQ

We ran

    ~/lcm/scripts/auto_barcode/barcode_split_trim.pl --barcode ../BCfile2.txt --list ../AdaptersRemoved.fq





