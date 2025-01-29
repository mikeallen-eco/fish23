# A metabarcoding bioinformatics tutorial using obitools, Swarm, and CRABS

NOTE: this tutorial uses an older version of CRABS (v0.8) and is therefore a little out of date. I plan to archive this and update it to work with CRABS v1.0.7 soon. It also uses obitools v2, which can be hard to install as it requires an old version of Python. 

This is an example workflow, using the Rutgers Amarel cluster to analyze eDNA samples taken to characterize freshwater fish communities within New Jersey, USA. The molecular analysis amplified the 'MiFish' section of the 12S mitochondrial region using Illumina next-generation sequencing (metabarcoding). The workflow follows the obitools pipeline (Boyer et al., 2016) with some modifications. It is largely based on the sources listed below and also code generously shared by others (especially B. Kwait, O. Stringham, O. Wangensteen, D. Marquina). But if you find it helpful to cite the full compiled workflow, you can also cite this page as:

Allen, M. C. 2024. A metabarcoding bioinformatics tutorial using obitools, Swarm, and CRABS. URL: https://github.com/mikeallen-eco/fish23

<ins>Sources</ins>
1. Boyer, F., Mercier, C., Bonin, A., Le Bras, Y., Taberlet, P., & Coissac, E. (2016). obitools: A unix‐inspired software package for DNA metabarcoding. Molecular ecology resources, 16(1), 176-182.
2. the obitools "wolf tutorial": 
https://pythonhosted.org/OBITools/wolves.html
3. a helpful pipeline published online by D. Marquina:
https://metagusano.github.io/publications/Bioinformatic%20Pipeline%20For%20Metabarcoding.pdf
4. Swarm output processing code by O. Wangensteen:
https://github.com/metabarpark/R_scripts_metabarpark/blob/master/owi_recount_swarm
5. Gold, Z., Curd, E. E., Goodwin, K. D., Choi, E. S., Frable, B. W., Thompson, A. R., ... & Barber, P. H. (2021). Improving metabarcoding taxonomic assignment: A case study of fishes in a large marine ecosystem. Molecular ecology resources, 21(7), 2546-2564.
6. Jeunen, G. J., Dowle, E., Edgecombe, J., von Ammon, U., Gemmell, N. J., & Cross, H. (2023). CRABS—a software program to generate curated reference databases for metabarcoding sequencing data. Molecular Ecology Resources, 23(3), 725-738.
7. Mahé, F., Rognes, T., Quince, C., de Vargas, C., & Dunthorn, M. (2014). Swarm: robust and fast clustering method for amplicon-based studies. PeerJ, 2, e593.

# Connect to the cluster
1. Connect to the Amarel cluster  
Visit the user guide to learn how to get an account on Amarel etc. https://sites.google.com/view/cluster-user-guide?pli=1

Open the shell in OnDemand (https://ondemand.hpc.rutgers.edu), the terminal if on a Mac, or MobaXterm (or Putty or similar). Log into the cluster by typing the code below. Be sure to connect to the VPN first if you are off campus (using AnyConnect software) - otherwise, it won't work.
```
ssh YOURNETID@amarel.rutgers.edu
```
2. Start an interactive session so you don't accidentally use the login nodes for processing.
```
srun -p main -N 1 -c 20 -n 1 -t 05:00:00 --mem 100GB --pty /bin/bash  # type that into a terminal/console
```

# Install software
1. Download miniconda if you don't have it installed already (for managing software versions)
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda -V # check that miniconda is installed
```
2. create a new conda environment that runs python v 2.7 and activate it
You'll need to activate it each time you log in in order to use any obitools commands.
```
conda create --name obi2 python=2.7
conda activate obi2
```
3. install obitools, ecopcr, cutadapt, prinseq, seqtk, swarm, and ecoprimers
install obitools version 2, ecopcr (also requires python 2.7)
```
conda install -c bioconda obitools
conda install -c bioconda ecopcr
conda install bioconda::cutadapt
conda install prinseq
conda install seqtk
```
4. Make a conda environment with python 3.8 and install CRABS & dependencies. Instructions and more info here: https://github.com/gjeunen/reference_database_creator
```
conda create -n py38 python=3.8
conda activate py38
conda install -c bioconda cutadapt=4.4
conda install -c bioconda vsearch
conda install -c bioconda muscle
pip install argparse
conda install -c conda-forge biopython # (this automatically installed numpy 1.2 also)
conda install -c conda-forge tqdm
conda install -c conda-forge matplotlib
conda install -c anaconda pandas
# install crabs by cloning the github to some location on your directory outside of your project
git clone https://github.com/gjeunen/reference_database_creator.git
conda install -c bioconda ecoprimers
conda install -c bioconda swarm # note: might be better to do this in the python 3.8 environment
```
You'll need to run this next line each time you start CRABS or else add it to your bashhc file (google that) to tell bash where to open the crabs program. (Note: change the file path to the actual one to the reference_database_creator folder.)
```
export PATH="/file/path/to/reference_database_creator:$PATH"
```
5. Test that it worked
```
crabs -h
```
6. Install NCBI's BLAST+ suite of programs
These two links are helpful: https://www.ncbi.nlm.nih.gov/books/NBK569856/ https://www.ncbi.nlm.nih.gov/books/NBK52640/
```
# make directory outside of the main project folder called blast
mkdir blast
cd blast

# download the tar.gz file
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz

# unzip it
tar -xvf ncbi-blast-2.14.0+-x64-linux.tar.gz

# add this next part to the .bashhc file or run it each time you want to use blastn
export PATH="/file/path/to/blast/ncbi-blast-2.14.0+/bin:$PATH"
```

# Process Illumina sequences
### Demultiplex and download the sequencing data
1. Demultiplexing means to split up the reads into individual files by sample. You can do this with bioninformatics software on the cluster, but it is currently pretty easy to do right on the Princeton HTSEQ website. To demultiplex on the website, log in and navigate to your "assay" (sequencing run) page. Click the green "Create Custom Barcode File" button & upload an index file as a tab delimited file with 3 columns: a unique ID, index1 sequences, index2 sequences. The unique ID is a concatenation of the sample pool (e.g., "X"), and the fields Index_Plate and Index_Plate_Well from the file: nextera-dna-udi-lrm-samplesheet-iSeq-MiSeq-nextera-xt-set-a-b-c-d-2x151-384-samples.csv. I'm not sure where that file came from, but I think it is a standard list of Illumina index "barcodes". Paste the 3 columns into a spreadsheet software and save as a tab delimited file. The unique ID could be anything, as long as allows you to match up the sequencing samples with a sample info data file.  
2. Once you have uploaded the index file, go back to the assay page and click "Demultiplex Assay Data" and "enqueue" the demultiplexing job with pre-set defaults (1 mismatch allowed, default file output format).
3. You'll receive an email telling you when it is done. To download, follow the link that they email you and click on the green "Download fastq" button. There, you can set the "Format" check box to "Wget" and copy the code they provide. It should look something like this: 
```
wget -r -nH --cut-dirs=2 https://htseq.princeton.edu/tmp/SOMECODE/
```
4. On the cluster, navigate to a directory with at least 1TB of storage space to spare. Create a main directory which I'll call "fish23" and a subdirectory called "rawdata".
```
mkdir fish23
mkdir fish23/rawdata
cd fish23/rawdata
```
5. Now that you are in the rawdata folder, run the wget commands for each sequencing pool in the command line. The downloads take only about 30 s per pool. In between downloads for each pool, move the output text files into additional subdirectories so they don't get over-written. OK, now all of the sequencing data is safely on the cluster.

### Merge paired reads

6. In the rawdata folder, move the unmatched reads to a new subfolder called 'unmatched'.
7. Activate the conda environment that contains obitools. For example:
```
conda activate obi2b
```
8. Go to the directory above rawdata and use the nano command to make a new bash script called 'f01_illuminapairedend.sh'.
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
9. Navigate back to the rawdata folder and run the script using the sbatch command. This will merge the paired reads (labeled read-1 and read-4) within the folder into a single file per read with the suffix .ali.fastq. Note: this can take a while, so can also be done in batches simultaneously by adding the sequencing pool number in front of the * above. E.g., if the prefix for sequencing pool X was 3011, then you'd make file f01_illuminapairedend.X.sh with the following bits of code modified:
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
10. Make a new directory within the main project directory called "aligned" and move the fastq files into subfolders within in by sequencing pool (e.g., folder X, Y, and Z).
### Trim off adapters and project codes
11. Use nano to make a script in the main directory called f02_cutadapt.sh and paste the text below into it. This will trim adapters from the sequences creating new files of trimmed sequences with the suffix *.ali.cut.fastq. You will need to make sure you are using the right adapter & project code sequences. (Viewing the sequences in a text file and searching for the primers & their reverse complements can be helpful with understanding this.) The code below also runs prinseq.sh to remove sequences w/ > 21 N bases. Change the file path to the prinseq-lite.pl script before you run this. Note that there are other ways to do this too if that software doesn't work. Run this script from within the directory your files are in.
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
18. Run seqtk.sh to turn fastq files into fasta files
```
for f in `ls -1 *.ali.cut.n21.fastq | sed 's/.ali.cut.n21.fastq//' `
do
seqtk seq -a ${f}.ali.cut.n21.fastq > ${f}.ali.cut.n21.fasta
done
```
### Remove primers, etc.
19. This is done using the obitools command ngsfilter. Use nano to make bash scripts to run the following code for each sequencing pool: f03_ngsfilter.X.sh, f03_ngsfilter.Y.sh, and f03_ngsfilter.Z.sh. Change the sequence pool code (3011 in the example below) in the file as needed. Run the jobs using sbatch called from within the rawdata folder. This script takes about a minute per file in its current form.
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=ngsX
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100000
#SBATCH --time=0-5:00:00
#SBATCH -o %j_%N.out
#SBATCH -e %j_%N.err

for f in `ls -1 3011*.n21.fasta | sed 's/.ali.cut.n21.fasta//' | sed 's/3011__//' `
do
# make an index file for each sample (overwriting it each time)
echo fish23 ${f} aaaaa:ttttt TCTTGTCGGTAAAACTCGTGCCAGC CCATAGTGGGGTATCTAATCCCAGTTTG > ngs_index.X.txt

# add a fake tag onto each sequence so that ngsfilter can run
sed '/^>/ !{s/^/aaaaa/; s/$/aaaaa/}' 3011__${f}.ali.cut.n21.fasta > 3011__${f}.ali.cut.n21.tag.fasta

# run ngsfilter to remove primers/codes and add sample name to header
ngsfilter -t ngs_index.X.txt -u unidentified${f}.fasta 3011__${f}.ali.cut.n21.tag.fasta > 3011__${f}.ali.cut.n21.tag.ngs.fasta
done
```
20. Merge all files into one.
```
cat *.ngs.fasta > merged.fasta
```
21. Dereplicate (takes about 30-40 min). Use nano to make a bash script (f04_mergeuni.sh) with the following code and run using sbatch.
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=obiuni
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80000
#SBATCH --time=0-01:00:00
#SBATCH --output=slurm.justr.%N.%j.out
#SBATCH --error=slurm.justr.%N.%j.err

# dereplicate
obiuniq -m sample merged.fasta > merged.uni.fasta

# remove unecessary headers
obiannotate -k count -k merged_sample merged.uni.fasta > $$ ; mv $$ merged.uni.fasta

# filter based on read length and minimum read count
# obigrep -l 140 -L 190 -p 'count>=0' merged.uni.fasta > merged.uni.c0.140.190.fasta
  # uncomment the above line if you want to test version with no min count threshold
obigrep -l 140 -L 190 -p 'count>=10' merged.uni.fasta > merged.uni.c10.140.190.fasta

```
Check the output.
```
obistat -c count merged.uni.c10.140.190.fasta |  sort -nk1 | head -20
```
22. ONLY DO THIS IF YOU DON't WANT TO USE SWARM TO CLUSTER. Run obiclean.sh to denoise (default settings? it is a simple form of clustering) - final file merged.uniq.c10.l140.l190.cln.fasta
```
#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=obiuni
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80000
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm.justr.%N.%j.out
#SBATCH --error=slurm.justr.%N.%j.err


obiclean -s merged_sample -r 0.05 -H \
  merged.uni.c10.140.190.fasta > merged.uni.c10.140.190.cln.fasta
```
# Clustering into MOTUs
Optional: clustering sequences into MOTUs with SWARM (if doing this, don't run obiclean in # 22 above)
Move merged.uni.c10.140.190.fasta into a new folder called "cluster"
1. shorten the uniq IDs (use a one-character prefix for it to work with the tallying script below)
```
obiannotate --seq-rank merged.uni.c10.140.190.fasta | obiannotate --set-identifier '"F_%d"%seq_rank' > merged.uni.c10.140.190.sht.fasta

obitab -o merged.uni.c10.140.190.sht.fasta > merged.uni.c10.140.190.sht.tab
```
2. Convert into vsearch format.
```
obiannotate -k count merged.uni.c10.140.190.sht.fasta > merged.uni.c10.140.190.sht.vsc.fasta

sed -i -e 's/[[:space:]]count/;size/g' merged.uni.c10.140.190.sht.vsc.fasta
sed -i 's/;\+$//' merged.uni.c10.140.190.sht.vsc.fasta
```
3. Remove the chimeras
Make a script file (f05_chimera.sh) and run it using sbatch. This removes reads judged to be chimeras by vsearch algorithm.
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=chimera
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80000
#SBATCH --time=0-00:30:00
#SBATCH --output=slurm.justr.%N.%j.out
#SBATCH --error=slurm.justr.%N.%j.err

vsearch -sortbysize merged.uni.c10.140.190.sht.vsc.fasta -output merged.uni.c10.140.190.sht.vsc.srt.fasta

vsearch --uchime_denovo merged.uni.c10.140.190.sht.vsc.srt.fasta --sizeout --nonchimeras merged.uni.c10.140.190.sht.vsc.srt.chi.fasta --chimeras merged.uni.c10.140.190.sht.vsc.srt.chimeras.fasta --uchimeout uchimeout.txt

# make it into a single-line fasta file
awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' merged.uni.c10.140.190.sht.vsc.srt.chi.fasta > merged.uni.c10.140.190.sht.vsc.srt.chi.sin.fasta
```
4. Use Swarm to cluster reads into MOTUs 
-d = max distance in nucleotides for 2 sequences to be same; Swarm github advises d=1 in most cases.
This one doesn't take too long, so it can be run in the command line in an interactive session.
```
swarm -d 1 -f -z -t 1 -o merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_output -s merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_stats -w merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fasta merged.uni.c10.140.190.sht.vsc.srt.chi.sin.fasta
```
5. Create R script to recount SWARM reads by OTU and sample.
Make a new file using nano or a similar program, paste the code below into it, and save it in the main directory of your project as f06_recount_swarm.R. Make sure "fileswarm" and "filetab" show the correct file path to the Swarm and obitab outputs, respectively. Note that this code also excludes singletons (clusters with total read count = 1 across samples), but you may have already cut out sequences with count < some minimum count threshold earlier (e.g., 10 in obigrep). The .tab file referenced was created in step 19b above. The script was modified from one found here: https://github.com/metabarpark/R_scripts_metabarpark/blob/master/owi_recount_swarm
```
## Script for recount abundances of a dataset clustered using Swarm to a tabulated csv file
## The script will read two arguments from the command line: the input file name (output file from swarm)
## and the abundances tab file from obitab.
## The output will be a ".csv" file with recalculated abundances per MOTU and per sample.
## By Owen S. Wangensteen - Project Metabarpark  2017
## Modified by Mike Allen, so that it would run with his data.

message("Please enter a cluster list output file obtained from SWARM and a tabulated counts file from obitab as fileswarm and filetab below")
  
  fileswarm <- "cluster/merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_output" #Must be an output file from swarm
  filetab <- "cluster/merged.uni.c10.140.190.sht.tab"
  outfile <-paste(fileswarm,".counts.csv",sep="")   
  

  
  min_reads <- 2
  
  get_swarm_size <- function(cadena="="){
    return(as.numeric(substr(cadena,gregexpr("=",cadena)[[1]][[1]]+1,nchar(cadena))))
  }
  
  # Read cluster list database
  message("Reading swarm database...")
  swarm_db <- readLines(fileswarm)
  total_swarms <- length(swarm_db)
  message("Cluster database read including ", total_swarms," total clusters.")
  
  # Calculate reads in each cluster 
  message("Calculating number of reads in each cluster")

  # detect if a semicolon exists at end of each cluster string and remove it if so (there never is?)
  clusters <- strsplit(swarm_db," ")
  for (i in 1:total_swarms) for (j in 1:length(clusters[[i]])) if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))==";"){
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-1)
  }
  
  # count total cluster size and make reduced swarm db removing clusters with less than threshold number of reads
  cluster_reads  <- NULL
  for (i in 1:total_swarms) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
  swarm_db_reduced <- swarm_db[cluster_reads>=min_reads]

  # make reduced cluster db into a list
  clusters <- strsplit(swarm_db_reduced," ")
  total_swarms_reduced <- length(swarm_db_reduced)
  
  # get cluster IDs from reduced swarm database
    # note: edited to autodetect number of sequence ID characters (because obitools seq_rank makes variable length IDs)
  id <- NULL
  for (i in 1:total_swarms_reduced) for (j in 1:length(clusters[[i]])) {
    num_char_id <- nchar(strsplit(clusters[[i]][[j]], split = ";")[[1]][[1]])
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,num_char_id)  
    id[i] <- clusters[[i]][1]
  }
  
  names(clusters) <- id
  
  message("Kept only ", total_swarms_reduced," clusters of size greater than or equal to ",min_reads," reads.")
  necesarios <- unlist(clusters, use.names=F)
  
  # Read counts database and keep only the needed clusters
  message("Reading tabulated database. This could take a while...")
  db <- read.table(filetab,sep="\t",head=T)
  numseqs <- nrow(db)
  db <- db[db$id %in% necesarios,]
  numseqs_reduced <- nrow(db)
  samples <- length(names(db)[substr(names(db),1,6)=="sample"])
  message("Database read including ", numseqs," total different sequences and ",samples," samples.")
  message("Kept only ", numseqs_reduced," sequences for calculations.") 
  # note: can be less if you removed rare clusters above and due to chimera removal
  
  db.total <- merge(data.frame(id),db,by="id") #Con esto se queda solo con los heads
  id <- db.total$id
  numclust <- nrow(db.total)
  
  for (fila in 1:numclust){
    head <- id[fila]
    tails <- unlist(clusters[names(clusters)==head])
    db.reduced <- db[db$id %in% tails,]
    suma <- colSums(db.reduced[,substr(names(db.total),1,6)=="sample"])
    db.total[fila,substr(names(db.total),1,6)=="sample"] <- suma
    db.total$cluster_weight[fila] <- nrow(db.reduced)
    message("Cluster ", fila, " / ",numclust, " ready, including ", db.total$cluster_weight[fila]," sequences.","\r",appendLF = FALSE)
  }
  db.total$total_reads <- rowSums(db.total[,substr(names(db.total),1,6)=="sample"])
  names(db.total[substr(names(db.total),1,6)=="sample"]) <- substr(names(db.total[substr(names(db.total),1,6)=="sample"]),8,nchar(names(db.total[substr(names(db.total),1,6)=="sample"])))
  write.table(db.total[,c(1:(ncol(db.total)-3),(ncol(db.total)-1):ncol(db.total),(ncol(db.total)-2))],outfile,sep=",",quote=F,row.names=F)
  message("File ", outfile, " written")
```
Change your conda environment to one that has R installed, run the Rscript, and change back to the conda environment with obitools installed.
```
conda activate ren423
Rscript f06_recount_swarm.R
conda activate obi2b
```
6. Fix header of the Swarm output fasta file 
Add space between MOTU name and the size variable.
```
sed 's/;size=/ size=/g' merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fasta > merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.fasta

# remove singletons from the SWARM fasta file: 
  # this is only necessary if you didn't set a minimum copy threshold in obigrep previously, e.g., <10 reads)

# obigrep -p 'size > 1' merged.uni.c0.140.190.sht.vsc.srt.chi.sin.sw1.fix.fasta > merged.uni.c0.140.190.sht.vsc.srt.chi.sin.sw1.fix.no1.fasta
```

# Create a reference database
Here we'll create a 'curated' reference database of annotated MiFish amplicons (i.e., with species names attached) to compare against our sample sequences. This is done here using CRABS software (Jeunen et al., 2022). Most of the code below is modified from the associated GitHub page: https://github.com/gjeunen/reference_database_creator

### Import and process private reference sequence data
Get CUSTOM reference sequences into crabs from your own fasta file (e.g., the Rees fish mitochondria bioproject)
Note: the first item in the fasta header has to be just an accession number or a species name (see crabs github). We will process these files all the way through taxonomy assignment first in the command line so we can see clearly what happens to them.

1. Make directories within the main one for all the raw reference database downloads.
```
mkdir fish
mkdir mam
mkdir herps
```
2. activate the conda environment with CRABS installed and import first fasta file of custom sequences into CRABS format.
Our first file is fish mitochondrial sequences from of Rees et al.. Note: this may be unnecessary as they all seem to have accession numbers from GenBank which we'll be accessing in the next section. Note the use of "--seq_header accession" as the sequence labels are accession numbers. In the next one, it will be species.
```
conda activate py38

crabs db_import --input rees.fasta --output custom1.fasta --seq_header accession --delim ' '
```
3. Use in silico pcr to get some amplicon (MiFish) sequences out of the fasta file you just imported. This will only extract amplicons that have primers attached, and only if the primers match the MiFish primers with error rate < 4.5.
```
crabs insilico_pcr \
--input custom1.fasta \
--output custom1.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC \
--rev CATAGTGGGGTATCTAATCCCAGTTTG \
--error 4.5
```
4. Use pga to exctract more sequences from the file
```
crabs pga --input custom1.fasta \
--output custom1.MiFish.pga60.fasta \
--database custom1.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG \
--speed slow --percid 0.60 --coverage 0.95 --filter_method strict
```
5. Import the second local database (step 1): remove "_12S" from the end of the fish sequence labels so that the label is an exact match for the species_name (with underscore in place of space).
```
sed 's/_12S//g' Fin_Clip_Species_List_NoFHCF.fasta > Fin_Clip_Species_List_NoFHCF.for.fasta 
```
6. Import the 2nd local sequence fasta into CRABS format. This one is from a Rutgers eDNA Lab 12S sequencing effort for fish and again has the species name (Genus_species) as the sequence label.
```
crabs db_import --input Fin_Clip_Species_List_NoFHCF.for.fasta --output custom2.fasta --seq_header species --delim ' '
```
7. Perform in silico pcr to get amplicon sequences out of the 2nd local fish sequence file. Note that the fwd primer is slightly different as the first part of the MiFish primer was cut off in the sequences file (i.e., this part was missing from the published MiFish forward primer: GTCGGTAAAAC).
```
crabs insilico_pcr \
--input custom2.fasta \
--output custom2.MiFish.fasta \
--fwd TCGTGCCAGC \
--rev CATAGTGGGGTATCTAATCCCAGTTTG \
--error 4.5
```
8. Use pga to exctract more sequences from the 2nd file. Again, we use the truncated MiFish forward primer.
```
crabs pga --input custom2.fasta \
--output custom2.MiFish.pga60.fasta \
--database custom2.MiFish.fasta \
--fwd TCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG \
--speed slow --percid 0.60 --coverage 0.95 --filter_method strict
```
9. Use sed to fix 2 typos in fish names in the 2nd local database
```
sed -i 's/Misgurnus_angullicaudatus/Misgurnus_anguillicaudatus/g' custom2.MiFish.pga60.fasta
sed -i 's/Anguillla_rostrata/Anguilla_rostrata/g' custom2.MiFish.pga60.fasta
```
10. This next part is just to test if the taxonomy assignment works. I.e., to make sure the sequences don't get put into the "missing" file and get lost when we do the final taxonomy assignment to the combined reference database.
```
crabs db_merge \
--output custom.MiFish.pga60.merged.test.fasta --uniq yes \
--input custom1.MiFish.pga60.fasta custom2.MiFish.pga60.fasta

crabs assign_tax \
--input custom.MiFish.pga60.merged.test.fasta  \
--output custom.MiFish.pga60.merged.test.tax.tsv \
--acc2tax nucl_gb.accession2taxid \
--taxid nodes.dmp \
--name names.dmp \
--missing missing_custom.MiFish.pga60.merged.test.tax.tsv

# remove the test files once they have been checked
rm *.test.*
```
### Get and process reference sequences from public databases
11. Download mitochondrial sequences from GenBank, MitoFish, etc. for fish, mammals, and 'herps' (birds, amphibians, reptiles) separately. First create an .sh script called ref01_mifish_fish.sh in your main directory and paste the following into it. Note that you'll need to either customize the code related to the custom reference databases below or else comment them out. GenBank search term tips to customize it (e.g., could just download relevant families): https://otagomohio.github.io/hacky2021/sessions/1005_ncbi/. Note: you'll need to change "your@email.edu" to your own email address to download from GenBank. Start with the fish:
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=MiFish
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=95GB
#SBATCH --time=1-10:00:00
#SBATCH -o %N_%j.mifish.out
#SBATCH -e %N_%j.mifish.err

# download all mitofish sequences
crabs db_download \
--source mitofish \
--output CRABS_Fish_mitofish_dl.fasta \
--keep_original no

# download all Actinopterygii mito sequences from NCBI
crabs db_download \
--source ncbi \
--database nucleotide \
--query '("Cladistia"[Organism] OR Cladistia[All Fields] OR "Teleostei"[Organism] OR Teleostei[All Fields] OR "Chondrostei"[Organism] OR Chondrostei[All Fields] OR "Holostei"[Organism] OR Holostei[All Fields]) AND mitochondrion[filter]' \
--output CRABS_FishActinopterygii_ncbi_dl.fasta \
--keep_original no \
--email your@email.edu \
--batchsize 5000

# download all non-Actinopterygii fish mito sequences from NCBI
crabs db_download \
--source ncbi \
--database nucleotide \
--query '("Cyclostomata"[Organism] OR Cyclostomata[All Fields] OR "Gnathostomata"[Organism] OR Gnathostomata[All Fields]) AND mitochondrion[filter] NOT ("Actinopterygii"[Organism] OR "Tetrapoda"[Organism])' \
--output CRABS_FishOthers_ncbi_dl.fasta \
--keep_original no \
--email your@email.edu \
--batchsize 5000

# merge NCBI download files into one file (to simplify pga analysis below) - also add in custom local sequences
crabs db_merge \
--output CRABS_Fish_ncbi_dl.fasta --uniq yes \
--input CRABS_FishActinopterygii_ncbi_dl.fasta CRABS_FishOthers_ncbi_dl.fasta

# merge all fish files into one file to perform in silico pcr on
crabs db_merge \
--output CRABS_Fish_ncbimito_dl.fasta --uniq yes \
--input CRABS_Fish_ncbi_dl.fasta CRABS_Fish_mitofish_dl.fasta

# in silico pcr to get sequences out of combined fish sequence file
crabs insilico_pcr \
--input CRABS_Fish_ncbimito_dl.fasta \
--output CRABS_Fish_ncbimito_dl.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC \
--rev CATAGTGGGGTATCTAATCCCAGTTTG \
--error 4.5

# pga to exctract more sequences from the combined NCBI data
crabs pga --input CRABS_Fish_ncbi_dl.fasta \
--output CRABS_Fish_ncbi_dl.MiFish.pga60.fasta \
--database CRABS_Fish_ncbimito_dl.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG \
--speed slow --percid 0.60 --coverage 0.95 --filter_method strict

# pga to exctract more sequences from MitoFish data
crabs pga --input CRABS_Fish_mitofish_dl.fasta \
--output CRABS_Fish_mitofish_dl.MiFish.pga60.fasta \
--database CRABS_Fish_ncbimito_dl.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG \
--speed slow --percid 0.60 --coverage 0.95 --filter_method strict

# merge the two pga files together, creating final database of amplicons
crabs db_merge \
--output CRABS_Fish_ncbimito_dl.MiFish.pga60.fasta \
--uniq yes \
--input CRABS_Fish_mitofish_dl.MiFish.pga60.fasta \
CRABS_Fish_ncbi_dl.MiFish.pga60.fasta

```
12. Now make one for mammals (mifish_mam.sh):
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=mamMiFish
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=95GB
#SBATCH --time=1-10:00:00
#SBATCH -o %N_%j.out
#SBATCH -e %N_%j.err

# download shorter mt sequences containing 12S from all mammals (including humans)
crabs db_download \
--source ncbi \
--database nucleotide \
--query '("Prototheria"[Organism] OR Prototheria[All Fields] OR "Theria"[Organism] OR Theria[All Fields]) AND mitochondrion[filter] AND ("1"[SLEN] : "10000"[SLEN])' \
--output CRABS_Mammals_ncbi_dlA.fasta \
--keep_original no \
--email your@email.edu \
--batchsize 5000

# download longer mt sequences containing 12S from all mammals except humans (long-read human mt genomes = very large file)
crabs db_download \
--source ncbi \
--database nucleotide \
--query '("Prototheria"[Organism] OR Prototheria[All Fields] OR "Theria"[Organism] OR Theria[All Fields]) AND mitochondrion[filter] AND ("10000"[SLEN] : "22000"[SLEN])' \
--output CRABS_Mammals_ncbi_dlB.fasta \
--keep_original no \
--email your@email.edu \
--batchsize 5000

# merge the 2 raw sequence download files together
crabs db_merge \
--output CRABS_Mammals_ncbi_dl.fasta --uniq yes \
--input CRABS_Mammals_ncbi_dlA.fasta CRABS_Mammals_ncbi_dlB.fasta 

# extract mifish sequences using in silico pcr
crabs insilico_pcr \
--input CRABS_Mammals_ncbi_dl.fasta \
--output CRABS_Mammals_ncbi_dl.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC \
--rev CATAGTGGGGTATCTAATCCCAGTTTG \
--error 4.5

# extract more sequences from the raw sequence database using pga
crabs pga --input CRABS_Mammals_ncbi_dl.fasta \
--output CRABS_Mammals_ncbi_dl.MiFish.pga60.fasta \
--database CRABS_Mammals_ncbi_dl.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG \
--speed slow --percid 0.60 --coverage 0.95 --filter_method strict

```
13. Now make one for birds/fish/amphibians (mifish_herps.sh):
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=herpMiFish
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=95GB
#SBATCH --time=1-10:00:00
#SBATCH -o %N_%j.out
#SBATCH -e %N_%j.err

# download all non-Mammal Tetrapod mito sequences from NCBI
crabs db_download \
--source ncbi \
--database nucleotide \
--query '(("Tetrapoda"[Organism] OR Tetrapoda[All Fields]) AND mitochondrion[filter] NOT "Mammalia"[Organism]' \
--output CRABS_Herps_ncbi_dl.fasta \
--keep_original no \
--email your@email.edu \
--batchsize 5000

# pcr to get MiFish sequences out of sequence file
crabs insilico_pcr \
--input CRABS_Herps_ncbi_dl.fasta \
--output CRABS_Herps_ncbi_dl.MiFish.fasta \
--fwd GTCGGTAAAACTCGTGCCAGC \
--rev CATAGTGGGGTATCTAATCCCAGTTTG \
--error 4.5

# extract more sequences from data download files using pga
crabs pga --input CRABS_Herps_ncbi_dl.fasta --output CRABS_Herps_ncbi_dl.MiFish.pga60.fasta \
--database CRABS_Herps_ncbi_dl.MiFish.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG \
--speed slow --percid 0.60 --coverage 0.95 --filter_method strict

```
14. Use sbatch to run each of the scripts you just made from within their respective folder: fish, mam, or herps.
```
cd fish
sbatch ../ref01_mifish_fish.sh
cd herps
sbatch ../ref02_mifish_herps.sh
cd mam
sbatch ../ref03_mifish_mam.sh
```
15. Combine all reference databases and add taxonomy
```
conda activate py38

# make directory for taxonomy files (under main directory)
mkdir crabtax
cd crabtax

# download taxonomy files to it
crabs db_download --source taxonomy

# navigate back to main directory
cd ..
mkdir refdb

# merge NCBI/mitofish/local pga (amplicon) files for all taxa into one file for taxonomy assignment
  # doing this in 2 steps so that --uniq yes can be used for those with unique accession numbers
    # 2nd custom sequences file imported into CRABS by species name had potential for duplicate names (even if diff seq)
      # maybe this isn't necessary if CRABS makes a unique accession ID for species-based imports
crabs db_merge \
--output refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.step1.fasta --uniq yes \
--input herps/CRABS_Herps_ncbi_dl.MiFish.pga60.fasta mam/CRABS_Mammals_ncbi_dl.MiFish.pga60.fasta fish/CRABS_Fish_ncbimito_dl.MiFish.pga60.fasta fish/custom1.MiFish.pga60.fasta 

crabs db_merge \
--output refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.fasta \
--input refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.step1.fasta fish/custom2.MiFish.pga60.fasta

rm refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.step1.fasta

# create final amplicon table with taxonomy
# add records that didn't match a taxon to the "missing" file to be reviewed later 
crabs assign_tax \
--input refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.fasta  \
--output refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.tsv \
--acc2tax crabtax/nucl_gb.accession2taxid \
--taxid crabtax/nodes.dmp \
--name crabtax/names.dmp \
--missing refdb/missing_CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.tsv

# dereplicate
crabs dereplicate --input refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.tsv \
--output refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.tsv --method uniq_species

# cleanup
crabs seq_cleanup --input refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.tsv \
--output refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.tsv \
--minlen 100 --maxlen 300 --maxns 10 --enviro yes --species yes

```
###Subset to include only species expected to occur in the area
The approach of taxonomic assignment based on a geographically localized reference database was adopted based on the work of Gold et al., (2021). It prevents many "genus-level" taxonomic assignments that can incorrectly result from lowest-common ancestor algorithms (e.g., ecotag) when the reference database includes species that are out of range; specifically, when there are out-of-range species that are closely related evolutionarily and therefore share very similar mitochondrial sequences (e.g., the MiFish 12S sequence). You can always assess later whether you think an out-of-range species is a plausible ID given the data by performing a blast search against a larger database later (see that section below). 

16. Make a version that is a subset of only local NJ species. First, make a text file with a list of all the vertebrate species found in New Jersey. This list was compiled from 4 websites: 
https://dep.nj.gov/njfw/fishing/freshwater/freshwater-fish-of-new-jersey/
https://dep.nj.gov/njfw/education/online-field-guide-for-reptiles-and-amphibians/
https://www.nj.gov/dep/fgw/chkbirds.htm
https://www.nj.gov/dep/fgw/chkmamls.htm

The lists were pasted into a csv file and processed using the make_nj_vert_list.R script to convert into NCBI taxonomy. This uses an 'align_taxonomy' function co-written with O. Stringham. Both scripts are found in the "R" directory of this repository. Paste the final list of local target species in NCBI format into a new text file created on the cluster using nano. Here we'll make one called nj_vertebrates_ncbi.txt with all the vertebrates found in New Jersey (except saltwater fish species):
```
Glyptemys muhlenbergii
Graptemys geographica
Sternotherus odoratus
Chelydra serpentina
Malaclemys terrapin
Terrapene carolina
Kinosternon subrubrum
Chrysemys picta
Apalone spinifera
Pseudemys rubriventris
Trachemys scripta
Clemmys guttata
Glyptemys insculpta
Plestiodon fasciatus
Scincella lateralis
Sceloporus undulatus
Pantherophis obsoletus
Pantherophis guttatus
Thamnophis sirtalis
Heterodon platirhinos
Lampropeltis getula
Lampropeltis triangulum
Thamnophis saurita
Virginia valeriae
Carphophis amoenus
Coluber constrictor
Storeria dekayi
Agkistrodon contortrix
Pituophis melanoleucus
Storeria occipitomaculata
Diadophis punctatus
Cemophora coccinea
Nerodia sipedon
Regina septemvittata
Opheodrys aestivus
Opheodrys vernalis
Crotalus horridus
Ambystoma laterale
Pseudotriton montanus
Ambystoma tigrinum
Hemidactylium scutatum
Ambystoma jeffersonianum
Eurycea longicauda
Ambystoma opacum
Desmognathus ochrophaeus
Desmognathus fuscus
Pseudotriton ruber
Plethodon glutinosus
Gyrinophilus porphyriticus
Eurycea bislineata
Plethodon cinereus
Notophthalmus viridescens
Ambystoma maculatum
Anaxyrus americanus
Aquarana catesbeiana
Lithobates virgatipes
Scaphiopus holbrookii
Anaxyrus woodhousii
Lithobates clamitans
Pseudacris triseriata
Acris crepitans
Dryophytes versicolor
Pseudacris crucifer
Lithobates palustris
Dryophytes andersonii
Dryophytes chrysoscelis
Lithobates sphenocephalus
Lithobates sylvaticus
Gavia stellata
Gavia immer
Podilymbus podiceps
Podiceps auritus
Podiceps grisegena
Calonectris diomedea
Ardenna gravis
Ardenna grisea
Oceanites oceanicus
Oceanodroma leucorhoa
Morus bassanus
Pelecanus occidentalis
Phalacrocorax carbo
Phalacrocorax auritus
Ixobrychus exilis
Ardea herodias
Ardea alba
Egretta thula
Egretta caerulea
Egretta tricolor
Bubulcus ibis
Butorides striata
Nycticorax nycticorax
Plegadis falcinellus
Dendrocygna bicolor
Cygnus columbianus
Cygnus olor
Anser caerulescens
Branta bernicla
Branta canadensis
Aix sponsa
Anas crecca
Anas rubripes
Anas platyrhynchos
Anas acuta
Spatula discors
Mareca strepera
Mareca penelope
Mareca americana
Aythya valisineria
Aythya americana
Aythya collaris
Aythya marila
Aythya affinis
Somateria mollissima
Somateria spectabilis
Histrionicus histrionicus
Clangula hyemalis
Melanitta nigra
Melanitta perspicillata
Melanitta fusca
Bucephala clangula
Bucephala albeola
Lophodytes cucullatus
Mergus merganser
Mergus serrator
Oxyura jamaicensis
Coragyps atratus
Cathartes aura
Pandion haliaetus
Ictinia mississippiensis
Haliaeetus leucocephalus
Circus cyaneus
Accipiter striatus
Accipiter cooperii
Accipiter gentilis
Buteo lineatus
Buteo platypterus
Buteo jamaicensis
Buteo lagopus
Aquila chrysaetos
Falco sparverius
Falco columbarius
Falco peregrinus
Phasianus colchicus
Bonasa umbellus
Meleagris gallopavo
Alectoris graeca
Colinus virginianus
Coturnicops noveboracensis
Laterallus jamaicensis
Rallus longirostris
Rallus elegans
Porzana carolina
Gallinula chloropus
Fulica americana
Pluvialis squatarola
Pluvialis dominica
Charadrius semipalmatus
Charadrius melodus
Charadrius vociferus
Haematopus palliatus
Himantopus mexicanus
Recurvirostra americana
Tringa melanoleuca
Tringa flavipes
Tringa solitaria
Tringa semipalmata
Actitis macularius
Bartramia longicauda
Numenius phaeopus
Limosa haemastica
Limosa fedoa
Arenaria interpres
Calidris canutus
Calidris alba
Calidris pusilla
Calidris mauri
Calidris minutilla
Calidris fuscicollis
Calidris bairdii
Calidris melanotos
Calidris maritima
Calidris alpina
Calidris ferruginea
Calidris himantopus
Calidris subruficollis
Calidris pugnax
Limnodromus griseus
Limnodromus scolopaceus
Gallinago gallinago
Phalaropus tricolor
Phalaropus lobatus
Stercorarius pomarinus
Stercorarius parasiticus
Stercorarius longicaudus
Leucophaeus atricilla
Hydrocoloeus minutus
Chroicocephalus ridibundus
Larus philadelphia
Larus delawarensis
Larus argentatus
Larus glaucoides
Larus fuscus
Larus hyperboreus
Larus marinus
Rissa tridactyla
Gelochelidon nilotica
Hydroprogne caspia
Thalasseus maximus
Sterna dougallii
Sterna hirundo
Sterna forsteri
Sternula antillarum
Rynchops niger
Alle alle
Uria lomvia
Alca torda
Columba livia
Zenaida macroura
Coccyzus erythropthalmus
Coccyzus americanus
Tyto alba
Megascops asio
Bubo virginianus
Bubo scandiacus
Strix varia
Asio otus
Asio flammeus
Aegolius acadicus
Chordeiles minor
Antrostomus carolinensis
Antrostomus vociferus
Chaetura pelagica
Archilochus colubris
Megaceryle alcyon
Melanerpes erythrocephalus
Melanerpes carolinus
Sphyrapicus varius
Dryobates pubescens
Picoides villosus
Colaptes auratus
Dryocopus pileatus
Contopus virens
Empidonax flaviventris
Empidonax virescens
Empidonax alnorum
Empidonax traillii
Empidonax minimus
Sayornis phoebe
Myiarchus crinitus
Tyrannus verticalis
Tyrannus tyrannus
Eremophila alpestris
Progne subis
Tachycineta bicolor
Stelgidopteryx serripennis
Riparia riparia
Petrochelidon pyrrhonota
Hirundo rustica
Cyanocitta cristata
Corvus brachyrhynchos
Corvus ossifragus
Corvus corax
Poecile carolinensis
Poecile hudsonicus
Baeolophus bicolor
Sitta canadensis
Sitta carolinensis
Certhia americana
Thryothorus ludovicianus
Troglogytes aedon
Troglodytes troglodytes
Cistothorus platensis
Cistothorus palustris
Regulus satrapa
Regulus calendula
Polioptila caerulea
Sialia sialis
Catharus fuscescens
Catharus minimus
Catharus ustulatus
Catharus guttatus
Hylocichla mustelina
Turdus migratorius
Dumetella carolinensis
Mimus polyglottos
Toxostoma rufum
Anthus rubescens
Bombycilla cedrorum
Lanius excubitor
Lanius ludovicianus
Sturnus vulgaris
Vireo griseus
Vireo solitarius
Vireo flavifrons
Vireo gilvus
Vireo philadelphicus
Vireo olivaceus
Vermivora cyanoptera
Vermivora chrysoptera
Leiothlypis peregrina
Leiothlypis celata
Setophaga americana
Setophaga petechia
Setophaga pensylvanica
Setophaga magnolia
Setophaga tigrina
Setophaga caerulescens
Setophaga coronata
Setophaga virens
Setophaga fusca
Setophaga dominica
Setophaga pinus
Setophaga discolor
Setophaga palmarum
Setophaga castanea
Setophaga striata
Setophaga cerulea
Mniotilta varia
Setophaga ruticilla
Protonotaria citrea
Helmitheros vermivorum
Seiurus aurocapilla
Parkesia noveboracensis
Parkesia motacilla
Geothlypis formosa
Oporornis agilis
Geothlypis philadelphia
Geothlypis trichas
Setophaga citrina
Cardellina pusilla
Cardellina canadensis
Icteria virens
Piranga rubra
Piranga olivacea
Cardinalis cardinalis
Pheucticus ludovicianus
Passerina caerulea
Passerina cyanea
Spiza americana
Pipilo erythrophthalmus
Spizelloides arborea
Spizella passerina
Spizella pusilla
Pooecetes gramineus
Chondestes grammacus
Passerculus sandwichensis
Ammodramus savannarum
Centronyx henslowii
Ammospiza caudacuta
Ammodramus nelsoni
Passerella iliaca
Melospiza melodia
Melospiza lincolnii
Melospiza georgiana
Zonotrichia albicollis
Zonotrichia leucophrys
Junco hyemalis
Calcarius lapponicus
Dolichonyx oryzivorus
Agelaius phoeniceus
Sturnella magna
Euphagus carolinus
Quiscalus major
Quiscalus quiscula
Molothrus ater
Icterus spurius
Icterus galbula
Pinicola enucleator
Haemorhous purpureus
Haemorhous mexicanus
Loxia curvirostra
Loxia leucoptera
Acanthis flammea
Spinus pinus
Spinus tristis
Passer domesticus
Alosa pseudoharengus
Lethenteron appendix
Anguilla rostrata
Alosa sapidissima
Monopterus albus
Salmo salar
Acipenser oxyrinchus
Fundulus diaphanus
Enneacanthus obesus
Hypophthalmichthys nobilis
Ameiurus melas
Pomoxis nigromaculatus
Enneacanthus chaetodon
Rhinichthys atratulus
Alosa aestivalis
Lepomis macrochirus
Enneacanthus gloriosus
Amia calva
Notropis bifrenatus
Culaea inconstans
Salvelinus fontinalis
Ameiurus nebulosus
Salmo trutta
Esox niger
Ictalurus punctatus
Notropis amoenus
Cyprinus carpio
Luxilus cornutus
Semotilus atromaculatus
Erimyzon oblongus
Exoglossum maxillingua
Gambusia holbrooki
Umbra pygmaea
Hybognathus regius
Semotilus corporalis
Pimephales promelas
Pylodictis olivaris
Dorosoma cepedianum
Notemigonus crysoleucas
Carassius auratus
Ctenopharyngodon idella
Lepomis cyanellus
Alosa mediocris
Trinectes maculatus
Notropis chalybaeus
Salvelinus namaycush
Micropterus salmoides
Rhinichthys cataractae
Lepisosteus osseus
Noturus insignis
Acantharchus pomotis
Fundulus heteroclitus
Esox masquinongy
Hypentelium nigricans
Esox lucius
Channa argus
Misgurnus anguillicaudatus
Aphredoderus sayanus
Lepomis gibbosus
Carpiodes cyprinus
Osmerus mordax
Oncorhynchus mykiss
Lepomis auritus
Esox americanus
Ambloplites rupestris
Cyprinella analostana
Petromyzon marinus
Percina peltata
Acipenser brevirostrum
Cottus cognatus
Micropterus dolomieu
Cyprinella spiloptera
Notropis hudsonius
Morone saxatilis
Notropis procne
Etheostoma fusiforme
Noturus gyrinus
Etheostoma olmstedi
Sander vitreus
Lepomis gulosus
Gambusia affinis
Ameiurus catus
Pomoxis annularis
Morone americana
Catostomus commersonii
Ameiurus natalis
Perca flavescens
Botaurus lentiginosus
Nyctanassa violacea
Spatula clypeata
Elanoides forficatus
Rallus limicola
Porphyrio martinica
Scolopax minor
Phalaropus fulicarius
Chlidonias niger
Contopus cooperi
Poecile atricapillus
Leiothlypis ruficapilla
Ammospiza maritima
Plectrophenax nivalis
Hesperiphona vespertina
Didelphis marsupialis
Sorex cinereus
Sorex palustris
Sorex fumeus
Sorex dispar
Blarina brevicauda
Parascalops breweri
Scalopus aquaticus
Condylura cristata
Myotis lucifugus
Myotis sodalis
Myotis septentrionalis
Myotis leibii
Lasionycteris noctivagans
Perimyotis subflavus
Eptesicus fuscus
Lasiurus borealis
Dasypterus intermedius
Aeorestes cinereus
Sylvilagus floridanus
Sylvilagus transitionalis
Lepus capensis
Lepus californicus
Lepus townsendii
Tamias striatus
Marmota monax
Sciurus carolinensis
Tamiasciurus hudsonicus
Glaucomys volans
Glaucomys sabrinus
Myocastor coypus
Oryzomys palustris
Peromyscus leucopus
Neotoma floridana
Clethrionomys gapperi
Microtus pennsylvanicus
Microtus pinetorum
Ondatra zibethicus
Synaptomys cooperi
Rattus rattus
Rattus norvegicus
Mus musculus
Napaeozapus insignis
Zapus hudsonius
Erethizon dorsatum
Vulpes vulpes
Urocyon cinereoargenteus
Ursus americanus
Procyon lotor
Mustela erminea
Mustela frenata
Neogale vison
Mephitis mephitis
Lontra canadensis
Lynx rufus
Odocoileus virginianus
Phoca vitulina
Halichoerus grypus
Cystophora cristata
Ziphius cavirostris
Mesoplodon densirostris
Mesoplodon europaeus
Mesoplodon mirus
Physeter catodon
Kogia breviceps
Kogia sima
Delphinapterus leucas
Stenella frontalis
Stenella plagiodon
Stenella coeruleoalba
Delphinus delphis
Tursiops truncatus
Orcinus orca
Grampus griseus
Globicephala melas
Phocoena phocoena
Balaenoptera physalus
Balaenoptera borealis
Balaenoptera acutorostrata
Balaenoptera musculus
Megaptera novaeangliae
Eubalaena glacialis
Homo sapiens
Bos taurus
Ovis aries
Gallus gallus
Capra hircus
Equus caballus
Sus scrofa
Sorex cinereus
Cryptotis parvus
Sorex hoyi
Castor canadensis
Canis latrans
Phoca groenlandica
Globicephala macrorhynchus
Canis lupus
Felis catus
Didelphis virginiana
```
17. Make a NJ subset of the reference database using the CRABS function db_subset and the list you just made. Note: same results if you swap the space within the binomials with an underscore.
```
crabs db_subset --input CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.tsv \
--output CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.njs.tsv \
--database nj_vertebrates_ncbi.txt --subset inclusion

# add header row to local NJ file
head CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.tsv -n 1 > tmp.tsv

cat tmp.tsv CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.njs.tsv > CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.njs2.tsv

rm tmp.tsv
```
18. Use the following R script to make the 2 final CRABS reference databases into fasta files: the full version and the NJ subset. The fasta will include a taxid for each entry that will work with obitools/ecotag taxonomy assignment. First, using nano, make an R script in the main directory called ref04_crabs2fasta.R with the following code in it:
```
# input the path to the full CRABS output reference file
crab_db_filepath_fulldb <- "refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.tsv"

# input the path to the NJ subset CRABS output reference file
crab_db_filepath_fulldb <- "refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.njs2.tsv"

# Define names for the full reference output file
output_ref_db_fulldb <- "refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.fasta"

# Define names for the NJ subset reference output file
output_ref_db <- "refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.njs.fasta"

# load packages
library(dplyr)
library(data.table)

###
# Make fasta file from full vertebrate reference db
###

# Load the database and add the taxonomic ranks
crab_db <- fread(crab_db_filepath_fulldb, sep="\t", head=T) %>%
  mutate(species = ifelse(species=="", "NA", species))

# Delete previous files
db_ref_file <- file(output_ref_db_fulldb,open = "wt")

closeAllConnections()

# Loop to write reference fasta file that will work with ecotag
showlines <- 1000

for (i in 1:nrow(crab_db)) {

        # write it in db_file
        output_file <- file(output_ref_db_fulldb, open = "at")
        writeLines(text = paste(">",as.character(crab_db$seqID[i])," species_name=",gsub(as.character(crab_db$species[i]),pattern="_",replacement = " "),
                                "; taxid=",as.integer(crab_db$taxid[i]),";",sep=""),
                   con = output_file)
        writeLines(text = gsub("-","",as.character(crab_db$sequence[i])), con = output_file)
        close(output_file)

  if (i %% showlines == 0) message(i,"/",nrow(crab_db)," sequences processed.","\r",appendLF = FALSE)
}

####
# Now make the NJ subset reference database fasta file
###

# Load the database and add the taxonomic ranks
crab_db <- fread(crab_db_filepath, sep="\t", head=T) %>%
  mutate(species = ifelse(species=="", "NA", species))

# Delete previous files
db_ref_file <- file(output_ref_db,open = "wt")

closeAllConnections()

# Loop to write reference fasta file that will work with ecotag
showlines <- 1000

for (i in 1:nrow(crab_db)) {

        # write it in db_file
        output_file <- file(output_ref_db, open = "at")
        writeLines(text = paste(">",as.character(crab_db$seqID[i])," species_name=",gsub(as.character(crab_db$species[i]),pattern="_",replacement = " "),
                                "; taxid=",as.integer(crab_db$taxid[i]),";",sep=""),
                   con = output_file)
        writeLines(text = gsub("-","",as.character(crab_db$sequence[i])), con = output_file)
        close(output_file)

  if (i %% showlines == 0) message(i,"/",nrow(crab_db)," sequences processed.","\r",appendLF = FALSE)
}
```
Next, edit the first two file paths to the correct ones for your reference database CRABS and the desired output file name. Now, run the script using the following commands from the main directory:
```
# activate conda environment with R installed
conda activate ren423

# run the script
Rscript ref04_crabs2fasta.R
```
# Taxonomy assignment with ecotag
This section reformats the reference database in obitools/EcoPCR format and performs taxonomy assignment using the ecotag algorithm. It might best be placed all in one .sh file and run as a batch.

1. Make a new directory within the main directory called taxo and download NCBI files. These will be used for assigning taxonomy with obitools/ecotag. Note: it can't hurt to redownload when in doubt as errors are possible when using older taxonomy databases.
```
mkdir taxo
cd taxo
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
```
2. Make an obitools/ecoPCR taxonomy database by running the following from within taxo directory:
```
obitaxonomy -t taxo -d taxo 
```
3. Assign taxonomy. Make a bash script in the main directory called f07_ecotag.sh and paste the code below into it. It will use obigrep to clean the reference database (keeping only those ID'd to species, genus, family level); obiuni to dereplicate reference sequences; and obiannotate to add unique IDs to the reference database. Somewhere in one of those commands, the full taxonomy is added to the header (including a taxid for each level). Lastly, the code assign taxonomy to the sequences in your samples using the ecotag "lowest common ancestor" algorithm. Navigate to the cluster" directory (the one with your sample sequences) and run it using sbatch from there. 
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=ecotag
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=95GB
#SBATCH --time=1-10:00:00
#SBATCH -o %N_%j.ecotag.out
#SBATCH -e %N_%j.ecotag.err

# add full taxonomy to reference header (?) and require species, genus, and family names
obigrep -d taxo/taxo --require-rank=species \
  --require-rank=genus --require-rank=family ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln.njs.fasta > ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln2.njs.fasta

# Dereplicate sequences (also adds full taxonomy to header, or maybe in previous step). 
obiuniq -d taxo/taxo \
  ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni.cln2.njs.fasta > ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni2.cln2.njs.fasta

# One more obigrep to ensure all dereplicated sequences have a taxid at the family level
obigrep -d taxo/taxo --require-rank=family \
  ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni2.cln2.njs.fasta > ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni2.cln3.njs.fasta

# Annotate reference database with unique IDs
obiannotate --uniq-id ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni2.cln3.fasta > ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni2.cln3.njs.ann.fasta

# assign taxonomy to main file using ecotag LCA algorithm
ecotag -d ../taxo/taxo -R ../refdb/CRABS_Vertebrates_ncbiMitofishLocal.MiFish.pga60.uni2.cln3.njs.ann.fasta \
merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.fasta > \
merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.fasta

# Remove unneeded fields from the header
    # note: I deleted extra stuff in the below script; see wolf tutorial for defaults
    # note: If running SWARM way, delete count as "size" is the important variable; otherwise DON'T delete count!
    
obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount \
--delete-tag=obiclean_count --delete-tag=obiclean_singletoncount \
--delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount \
--delete-tag=obiclean_head --delete-tag=taxid_by_db --delete-tag=obiclean_headcount \
--delete-tag=id_status --delete-tag=rank_by_db --delete-tag=count \
--delete-tag=order merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.fasta > \
merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.ann.fasta

# Sort file by total read count. 
# Note: if you skipped the steps for clustering with Swarm, then the first argument would be count instead of size

obisort -k size -r merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.ann.fasta >  \
merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.ann.srt.fasta

# Export MOTU read count table by running the following code in the command line.

obitab -o merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.ann.srt.fasta > \
FishIBI.2023.final.MOTU.table.MiFish.tsv
```
4. download the tsv file and use process_tabs.Rmd script to do post processing (checking negatives, etc). If you ran the SWARM pathway, then you'll also need to download the *count.csv file and join it to the tsv using the MOTU id (or the sequence). 

# Perform BLAST on sequences
Even with a well-made reference database, there is always the possibility that some valid sequences were left out, either excluded during the download and filtering steps, or during the geographic filtering step. Therefore, I find it is a good idea to blast all MOTUs against the full NCBI Eukaryotic database and manually examine the results. This rarely produces any big surprises, but it often reveals a few better matches to the same taxa, or reveals likely contaminant species (or maybe recent introductions!) present at low total read counts. For example, it has revealed ocean fish species in freshwater river samples, all with very few reads (< 30 total).

1. download the entire eukaryotic component of the blast database or nt_euk. It comes from here https://ftp.ncbi.nlm.nih.gov/blast/db/. I like to put it in a subfolder of the blast directory called "blastdb" created when installing blast. Note: the above only worked for me when I was in my obitools conda environment. Navigate to blastdb directory and run the following code as an sh file using sbatch: sbatch ../../../fish23/f08_dl_blast.sh
Not sure what I did to install perl etc. but perhaps using conda within my obitools environment.


```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=blastdl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=95GB
#SBATCH --time=1-10:00:00
#SBATCH -o %N_%j.blast.out
#SBATCH -e %N_%j.blast.err

perl ../bin/update_blastdb.pl --passive --decompress nt_euk
# perl ../bin/update_blastdb.pl --passive --decompress nt_prok

```
3. Blast the unassigned reads. Make a sh script called f09_blast.sh in the main directory with the following code in it. Navigate to the blastdb folder you made (where you downloaded all the blast databases) and run it from there. You need to be set to the directory of the blast db AND the taxonomy db when you run the bash script, otherwise, species will be NA. Found that info here: https://www.biostars.org/p/76551/. Run the script using something like this:
```
sbatch ../../../fish23/f09_blast.sh
```
Here is the contents of f09_blast.sh (you'll need to change file paths):
```
#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=blast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=95GB
#SBATCH --time=1-10:00:00
#SBATCH -o %N_%j.blast.out
#SBATCH -e %N_%j.blast.err

export PATH="/file/path/to/blast/ncbi-blast-2.14.0+/bin:$PATH"

# blast the eukaryotic database
blastn -query "/file/path/to/maindir/cluster/merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.ann.srt.fasta" \
-db /file/path/to/blast/ncbi-blast-2.14.0+/blastdb/nt_euk \
-num_threads 24 \
-outfmt "6 delim=, std qlen slen staxids sscinames scomnames sskingdoms" -max_target_seqs 50 \
-out "/file/path/to/maindir/refdb/all.MOTUs.blasted.txt"

# blast the prokaryotic database
# blastn -query "/file/path/to/maindir/cluster/merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fix.tax.njs.ann.srt.fasta" \
# -db /file/path/to/blast/ncbi-blast-2.14.0+/blastdb/nt_prok \
# -num_threads 24 \
# -outfmt "6 delim=, std qlen slen staxids sscinames scomnames sskingdoms" -max_target_seqs 50 \
# -out "/file/path/to/maindir/refdb/all.MOTUs.blasted.prok.txt"

```
# Process final output
Download the following files to your computer from the cluster and add them to the data folder of this repository: 
1. merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_output.counts.csv
2. FishIBI.2023.final.MOTU.table.MiFish.tsv
3. all.MOTUs.blasted.txt
4. all.MOTUs.blasted.prok.txt (if applicable; I didn't end up using this one)

Open the process_tabs.Rmd in this repository in RStudio and follow instructions in there to process the data the rest of the way.

    
