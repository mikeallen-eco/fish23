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
You'll need to activate it each time you log in in order to use any obitools commands.
```
conda create --name obi2 python=2.7
conda activate obi2
```
5. install obitools, ecopcr, cutadapt, prinseq, seqtk, swarm, and ecoprimers
install obitools version 2, ecopcr (also requires python 2.7)
```
conda install -c bioconda obitools
conda install -c bioconda ecopcr
conda install bioconda::cutadapt
conda install prinseq
conda install seqtk
```

Install swarm (I think install this one on a different conda env with python 3+)
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
9. On the cluster, navigate to a directory with at least 1TB of storage space to spare. Create a main directory which I'll call "fish23" and a subdirectory called "rawdata".
```
mkdir fish23
mkdir fish23/rawdata
cd fish23/rawdata
```
10. Now that you are in the rawdata folder, run the wget commands for each sequencing pool in the command line. The downloads take only about 30 s per pool. In between downloads for each pool, move the output text files into additional subdirectories so they don't get over-written. OK, now all of the sequencing data is safely on the cluster.

# Merge paired reads

11. In the rawdata folder, move the unmatched reads to a new subfolder called 'unmatched'.
12. Activate the conda environment that contains obitools. For example:
```
conda activate obi2b
```
13. Go to the directory above rawdata and use the nano command to make a new bash script called 'f01_illuminapairedend.sh'.
```
nano f01_illuminapairedend.sh
```
Paste the following inside the script file:
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=illuminapaired
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100000
#SBATCH --time=0-8:00:00
#SBATCH -o %j_%N.out
#SBATCH -e %j_%N.err

# loop to merge reads
for f in `ls -1 *-read-[1,4].fastq.gz | sed 's/-read-[1,4].fastq.gz//' `
do
illuminapairedend --score-min=40 -r ${f}-read-4.fastq.gz ${f}-read-1.fastq.gz > ${f}.fastq
done

# loop to keep only successfully aligned reads
for f in `ls -1 *.fastq | sed 's/.fastq//' `
do
obigrep -p 'mode!="joined"' ${f}.fastq > ${f}.ali.fastq
done
```
14. Navigate back to the rawdata folder and run the script using the sbatch command. This will merge the paired reads (labeled read-1 and read-4) within the folder into a single file per read with the suffix .ali.fastq. Note: this can take a while, so can also be done in batches simultaneously by adding the sequencing pool number in front of the * above. E.g., if the prefix for sequencing pool X was 3011, then you'd make file f01_illuminapairedend.X.sh with the following bits of code modified:
```
# modified part of first loop
ls -1 3011*-read-[1,4].fastq.gz | sed 's/-read-[1,4].fastq.gz//'

# modified part of second loop
ls -1 3011*.fastq | sed 's/.fastq//'
```
This is how to run the script using sbatch (again, run this from within the rawdata folder):
```
sbatch ../f01_illuminapairedend.X.sh
```
15. Make a new directory within the main project directory called "aligned" and move the fastq files into subfolders within in by sequencing pool (e.g., folder X, Y, and Z).
# Trim adapters and primers from reads
16. Use nano to make a script in the main directory called f02_cutadapt.sh and paste the text below into it. This will trim primers and adapters from the sequences creating new files of trimmed sequences with the suffix *.ali.cut.fastq. It also runs prinseq.sh to remove sequences w/ > 21 N bases. Change the file path to the prinseq-lite.pl script before you run this. Run the script from within the directory your files are in.
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100000
#SBATCH --time=0-7:00:00
#SBATCH -o %j_%N.out
#SBATCH -e %j_%N.err

# loop to run cutadapt (remove any remaining adapters, trim low-quality reads)
for f in `ls -1 *.ali.fastq | sed 's/.ali.fastq//' `
do
cutadapt -q 30 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG  -g GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o ${f}.ali.cut.fastq ${f}.ali.fastq
done

# loop to remove sequences with more than 21 N bases
for f in `ls -1 *.ali.cut.fastq | sed 's/.ali.cut.fastq//' `
do
perl /YOURHOMEDIRECTORY/miniconda3/envs/obi2b/bin/prinseq-lite.pl -verbose -fastq ${f}.ali.cut.fastq -ns_max_n 21 -out_good ${f}.ali.cut.n21 -out_bad ${f}.ali.cut.y21n
done
```
17. Run fastqc for quality control
```
fastqc
```
The resulting html files should be downloaded and examinied to check that you actually got rid of adapters, low-qual bases, etc.
19. Run seqtk.sh to turn fastq files into fasta files
```
for f in `ls -1 *.ali.cut.fastq | sed 's/.ali.cut.fastq//' `
do
seqtk seq -a ${f}.ali.cut.fastq > ${f}.ali.cut.fasta
done
```
# Add sample name to sequence header

```
for f in `ls -1 *.n21.fasta | sed 's/.ali.cut.n21.fasta//' | sed 's/3011__//' `
do
obiannotate -S sample:${f} 3011__${f}.ali.cut.n21.fasta > ${f}.ali.cut.n21.ann.fasta
done
```

parts below here are placeholders that will be revised...
add fake indices (aaaaa:aaaaa) to sequences of each sample separately using indexS.sh scripts where S is sample name
    used custom R script to make each .sh file (make_bash.Rmd), uploaded to respective folders in "aligned" using OnDemand file explorer 
12. add sample name to each sequence header by running ngsS.sh scripts where S is sample name and indexS.txt is an input file with sample info
    used custom R script to make each .sh and .txt file (make_bash.Rmd)
# note: #11 & 12 might be avoidable with obiannotate: https://git.metabarcoding.org/obitools/obitools3/-/issues/127
    # also it may be possible to run ngsfilter with blank tags using --:-- but I forgot where I read that! 
13. make new directory: merged
    mkdir merged
14. copy the final sample files to the folder called merged
cp *.fasta ../merged # that would copy everything in your working directory ending in .fasta to merged
14a. run the concat.sh scripts to compile all samples into a single file; delete the individual sample file copies   
15. run obiuni.sh to collapse all reads to counts by sample
16. run obiannotate.sh to remove unnecessary header info
17. run obigrep2.sh, renaming file to merged.uniq.c10.l140.L190.fasta 
  17a. 	obigrep settings: >= 10 read count, for MiFish: length 140-190 (from Mjolnir documentation & MiFish paper)
18. ONLY DO THIS IF YOU DON't WANT TO USE SWARM TO CLUSTER. Run obiclean.sh to denoise (default settings? it is a simple form of clustering) - final file merged.uniq.c10.l140.l190.cln.fasta

