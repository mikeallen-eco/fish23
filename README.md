# fish23 (note: under construction)
fish metabarcoding bioinformatics workflow

# Connect to the cluster

1. Connect to the Amarel cluster  
Go to the shell in OnDemand, the terminal if on a Mac, or MobaXterm (or Putty or similar), log into the cluster and type the code below. Be sure to connect to the VPN first if you are off campus (using AnyConnect software) - otherwise, it won't work.
```
ssh YOURNETID@amarel.rutgers.edu
```
2. Start an interactive session so you don't accidentally use the login nodes for processing.
```
srun -p main -N 1 -c 2 -n 1 -t 05:00:00 --pty /bin/bash  # type that into a terminal/console
```

# Set up the required software
3. Download miniconda if you don't have it installed already (for managing software versions)
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda -V # check that miniconda is installed
```
4. create a new conda environment that runs python v 2.7 and activate it
You'll need to activate it each time you log in in order to use any obitools commands
```
conda create --name obi2 python=2.7
conda activate obi2
```
5. install obitools, ecopcr, swarm, and ecoprimers
install obitools version 2, ecopcr (also requires python 2.7)
```
conda install -c bioconda obitools
conda install -c bioconda ecopcr
```

install swarm (I think install this one on a different conda env with python 3+)
```
conda install -c bioconda swarm
```
install ecoprimers
```
conda install -c bioconda ecoprimers
```
# Demultiplex and download the sequencing data
6. Demultiplexing means to split up the reads into individual files by sample. You can do this with bioninformatics software on the cluster, but it is currently pretty easy to do right on the Princeton HTSEQ website. To demultiplex on the website, log in and navigate to your "assay" (sequencing run) page. Click the green "Create Custom Barcode File" button & upload an index file as a tab delimited file with 3 columns: a unique ID, index1 sequences, index2 sequences. The unique ID is a concatenation of the sample pool (e.g., "X"), and the fields Index_Plate and Index_Plate_Well from the file: nextera-dna-udi-lrm-samplesheet-iSeq-MiSeq-nextera-xt-set-a-b-c-d-2x151-384-samples.csv. I'm not sure where that file came from, but I think it is a standard list of Illumina index "barcodes". Paste the 3 columns into a spreadsheet software and save as a tab delimited file. The unique ID could be anything, as long as allows you to match up the sequencing samples with a sample info data file.  
7. Once you have uploaded the index file, go back to the assay page and click "Demultiplex Assay Data" and "enqueue" the demultiplexing job with pre-set defaults (1 mismatch allowed, default file output format).
8. You'll receive an email telling you when it is done. To download, follow the link that they email you and click on the green "Download fastq" button. There, you can set the "Format" check box to "Wget" and copy the code they provide. It should look something like this: 
```
wget -r -nH --cut-dirs=2 https://htseq.princeton.edu/tmp/SOMECODE/
```
9. On the cluster, navigate to a directory with at least 1TB to spare. Create a main directory which I'll call "fish23" and a subdirectory called "rawdata".
```
mkdir fish23
cd fish23
mkdir rawdata
cd rawdata
```
10. Now that you are in the raw data file, run the wget commands for each sequencing pool in the command line. The downloads take only about 30 s per pool. In between downloads for each pool, move the output text files into additional subdirectories so they don't get over-written. OK, now all of the sequencing data is safely on the cluster.

