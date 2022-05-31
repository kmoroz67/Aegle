

## Introduction
This pipeline is a tool for determining the clonotypes of T and B cell receptors,
as well as the main types of HLA class 1 and 2 alleles and gene expression in the sample according to RNA sequencing.

List of used bioinformatics tools:

- [MultiQC](https://github.com/ewels/MultiQC) - tool for RNA-seq data QC

- [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) - tool for converting SRA data to .fastq format;
is a faster version of the fastq-dump tool

- [MiXCR](https://github.com/milaboratory/mixcr) - versatile tool for fast and accurate analysis of T- and B cell repertoire based on NGS data.

- [Optitype](https://github.com/FRED-2/OptiType) - algorithm for HLA genotyping based on NGS data. Allows you to determine the alleles of HLA class 1.

- [Kallisto](https://chmi-sops.github.io/mydoc_kallisto.html) - Pseudoaligner for quick genes expression calculation. Returns expressions in TPM and counts.

- [Seq2HLA](https://github.com/TRON-Bioinformatics/seq2HLA) - In-silico method written in Python and R to determine HLA class 1/2 genotypes of a sample.

## Pipeline structure:  
Pipeline consists of several scripts:

- functions.py, which contains all the functions used in the pipeline
- geo_downloader.py, which contains the part of the pipeline for downloading and processing data from the GEO database (Gene Expression Omnibus)
- tgca_downloader.py, which contains the part of the pipeline for downloading and processing data from the GDC database (Genomic Data Commons)
- main.py, from which, depending on the parameter passed, either geo_downloader.py or tcga_downloader.py is called
- install_libs.sh, bash script for installation of required python libraries

## Set ups

**2. Set up fasterq-dump**  

Download archive (https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-centos_linux64.tar.gz) and unpack to the folder named "ncbi".  
	
	$ tar -xvzf sratoolkit.2.10.9-centos_linux64-cloud.tar.gz -C ./ncbi
	
Pull docker (if run at local machine/computer use sudo docker pull):
   	$ docker pull fred2/optitype
	
**3. Install dependencies**	

Run bash script to install libraries to your python virtual environment
	$ ./install_libs.sh
	
You also can install all dependencies using pip/pip3 in "hand mode"

**4. Install samtools**

	$ apt update
	$ apt install samtools

## Calculation set up
#### 1. Clone git repository

    $ git clone https://github.com/kmoroz67/Aegle.git

#### 2. Change directory to Aegle

    $ cd Aegle/

#### 3. Run pipeline  

Example of run from linux terminal:

!parameters in round bracers () is optional

    $ python main.py (-u <user_name>) --dbase tcga --n_start 0 --n_end 7

You can get a brief summary of the parameters used with the following command:

    $ python main.py -h

## Console interface
This tool is a console utility that supports the following options:

- -u / --user (optional) The name of the user in whose home folder the mixcr_calc folder will be created to store intermediate results.
By default, the name of the user on behalf of which the pipeline is launched is indicated.
	
- --dbase - the name of the database from which the data is downloaded and processed: geo or tcga
	
- --n_start - number of the first cpu from the range of cpus that will be used for MiXCR calculations

- --n_end - number of the last cpu from the range of cpus that will be used for MiXCR calculations

	The -h / --help flag is also available, displaying brief help about using the tool

## Results storage

Absolute path: **/<user>/calculated_datasets/...**

The downloaded data has the following structure:  dataset --> sample --> results folders

- hla1_results - folder with Optitype results 

- mixcr_results - folder with MiXCR results
	
- kallisto_results - folder with Kallisto results
	
- seq2hla_results - folder with Seq2HLA results

There is also a **tmp** folder that contains the following files: 

For data from the GEO portal:

- **ann_calc_slice.csv** - table of calculation queue   

	Run - run name   
	Sample - sample name  
	Dataset - dataset name
	
Layout of the data (paired or single) will be determined automatically.
	
- **download_process_samples.csv** - list of samples that were queued for processing by the pipeline

	Sample - sample name 
	
- **pipeline.log**
	
	Contains logs on the operation of the geo_downloader.py file
