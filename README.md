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
srun -p main -N 1 -c 20 -n 1 -t 05:00:00 --mem 100GB --pty /bin/bash  # type that into a terminal/console
```

# Set up the required software
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
You'll need to run this next line each time you start CRABS or else add it to your bashhc file (google that) to tell bash where to open the crabs program.
```
export PATH="/projects/f_deenr_1/mcallen/reference_database_creator:$PATH"
```
5. Test that it worked
```
crabs -h
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
### Trim adapters and primers from reads
11. Use nano to make a script in the main directory called f02_cutadapt.sh and paste the text below into it. This will trim primers and adapters from the sequences creating new files of trimmed sequences with the suffix *.ali.cut.fastq. It also runs prinseq.sh to remove sequences w/ > 21 N bases. Change the file path to the prinseq-lite.pl script before you run this. Run the script from within the directory your files are in.
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
Make a new file using nano or a similar program, paste the code below into it, and save it in the main directory of your project as f06_recount_swarm.R. Make sure "fileswarm" and "filetab" show the correct file path to the Swarm and obitab outputs, respectively. Note that this code also excludes singletons (clusters with total read count = 1 across samples), but you may have already cut out sequences with count < some minimum count threshold earlier (e.g., 10 in obigrep). The .tab file referenced was created in step 19b above.
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

# Create a custom reference database using CRABS

1. Make directories within the main one for all the raw reference database downloads.
```
mkdir fish
mkdir mam
mkdir herps
mkdir reftax
```
2. Curate any custom reference databases
Get CUSTOM reference sequences into crabs from your own fasta file (e.g., the Rees fish mitochondria bioproject)
Note: the first item in the fasta header has to be an accession number or a species name (see crabs github).

```
conda activate py38 # activate the conda environment with CRABS installed
crabs db_import --input rees.fasta --output custom1.fasta --seq_header accession --delim ' '
crabs db_import --input Fin_Clip_Species_List_NoFHCF.fasta --output custom2.fasta --seq_header species --delim ' '
# NOTE: to do: use sed to remove "_12S" from all the sequences so CRABS can read and assign taxonomy
```
Download mitochondrial sequences from GenBank, MitoFish, etc. for fish, mammals, and 'herps' (birds, amphibians, reptiles) separately. First create an .sh script called ref01_mifish_fish.sh in your main directory and paste the following into it. Note that you'll need to either customize the code related to the custom reference databases below or else comment them out. GenBank search term tips to customize it (e.g., could just download relevant families): https://otagomohio.github.io/hacky2021/sessions/1005_ncbi/. Note: you'll need to change "your@email.edu" to your own email address to download from GenBank. Start with the fish:

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
--input CRABS_FishActinopterygii_ncbi_dl.fasta CRABS_FishOthers_ncbi_dl.fasta custom1.fasta custom2.fasta

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

# merge the two pga files together, creating final database of 12SV5 amplicons
crabs db_merge \
--output CRABS_Fish_ncbimito_dl.MiFish.pga60.fasta \
--uniq yes \
--input CRABS_Fish_mitofish_dl.MiFish.pga60.fasta \
CRABS_Fish_ncbi_dl.MiFish.pga60.fasta

# download taxonomy files to it (comment out if already done)
crabs db_download --source taxonomy

# create final amplicon table with taxonomy
# add records that didn't match a taxon to the "missing" file to be reviewed later 
crabs assign_tax \
--input CRABS_Fish_ncbimito_dl.MiFish.pga60.fasta  \
--output CRABS_Fish_ncbimito_dl.MiFish.pga60.tax.tsv \
--acc2tax nucl_gb.accession2taxid \
--taxid nodes.dmp \
--name names.dmp \
--missing missing_CRABS_Fish_ncbimito_dl.MiFish.pga60.tax.tsv
```
Now make one for mammals (mifish_mam.sh):
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

# download taxonomy files to it (commented out since already done in fish folder)
# crabs db_download --source taxonomy

# create final amplicon table with taxonomy
# add records that didn't match a taxon to the "missing" file to be reviewed later 
crabs assign_tax \
--input CRABS_Mammals_ncbi_dl.MiFish.pga60.fasta \
--output CRABS_Mammals_ncbi_dl.MiFish.pga60.tax.tsv \
--acc2tax nucl_gb.accession2taxid \
--taxid nodes.dmp \
--name names.dmp \
--missing missing_CRABS_Mammals_ncbi_dl.MiFish.pga60.tax.tsv

```
Now make one for birds/fish/amphibians (mifish_herps.sh):
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

# download taxonomy files to it (commented out since already done in fish folder)
crabs db_download --source taxonomy

# create final amplicon table with taxonomy
# add records that didn't match a taxon to the "missing" file to be reviewed later 
crabs assign_tax \
--input CRABS_Herps_ncbi_dl.MiFish.pga60.fasta  \
--output CRABS_Herps_ncbi_dl.MiFish.pga60.tax.tsv \
--acc2tax nucl_gb.accession2taxid \
--taxid .nodes.dmp \
--name names.dmp \
--missing missing_CRABS_Herps_ncbi_dl.MiFish.pga60.tax.tsv

```
Use sbatch to run each of the scripts you just made from within their respective folder: fish, mam, or herps.
```
cd fish
sbatch ref01_mifish_fish.sh
```



# ... below here is really unfinished







22. assign taxonomy to the reference database
crabs assign_tax --input ncbi.bold.art.hum.mus.may23.pcr.pga.C.fasta --output ncbi.bold.art.hum.mus.may23.pcr.pga.C.tsv --acc2tax ../crabtaxo/nucl_gb.accession2taxid --taxid ../crabtaxo/nodes.dmp --name ../crabtaxo/names.dmp --missing missing_taxa.tsv
crabs assign_tax --input ncbi.bold.art.hum.mus.may23.pcr.pga.Z.fasta --output ncbi.bold.art.hum.mus.may23.pcr.pga.Z.tsv --acc2tax ../crabtaxo/nucl_gb.accession2taxid --taxid ../crabtaxo/nodes.dmp --name ../crabtaxo/names.dmp --missing missing_taxa.tsv
crabs assign_tax --input combined.fish.reference.may23.pcr.pga--output combined.fish.reference.may23.pcr.pga.tax.tsv --acc2tax ../crabtax/nucl_gb.accession2taxid --taxid ../crabtaxo/nodes.dmp --name ../crabtaxo/names.dmp --missing missing_taxa.tsv
23. remove extra duplicate columns added by assign_tax due to our splitting/recombining
    # note: counts the number of columns & keeps only the first 10 columns
    # note: might not be needed in all cases
awk '{print NF}' combined.fish.reference.may23.pcr.sin.pga.tax.tsv | sort -nu | tail -n 1
cut -f -10 combined.fish.reference.may23.pcr.sin.pga.tax.tsv > combined.fish.reference.may23.pcr.sin.pga.tax.fix.tsv
24. dereplicate reference sequences
crabs dereplicate --input ncbi.bold.art.hum.mus.may23.pcr.pga.tax.Z.tsv --output ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.Z.tsv --method uniq_species
crabs dereplicate --input combined.fish.reference.may23.pcr.sin.pga.tax.fix.tsv --output combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.tsv --method uniq_species
25. clean up remaining sequences
crabs seq_cleanup --input ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.C.tsv --output ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.fix.cln.C.tsv --minlen 25 --maxlen 500 --maxns 0 --enviro yes --species yes --nans 0
crabs seq_cleanup --input ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.Z.tsv --output ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln.Z.tsv --minlen 25 --maxlen 500 --maxns 0 --enviro yes --species yes --nans 0
crabs seq_cleanup --input ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.C.tsv --output ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln.C.tsv --minlen 25 --maxlen 500 --species yes
crabs seq_cleanup --input combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.tsv --output combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln.tsv --minlen 25 --maxlen 500 --species yes
26. use R to convert tsv to obitools-compatible fasta using script adapted from Owen W;
header should be ~: >GBMOR2153-19 species_name=Schistocerca americana; rank=species; origin = bld_gb; taxid=7009;
# the script "crabs2ecotag.Rmd" works for now, but only outside of the cluster
# NOTE: at this stage you can also filter your reference database based on geography (in R) see inspect_ref_db.R
    # for example, for fish, we extracted just those fish found in NY/NJ/PA + exotics + non-fish vertebrates
27. download taxonomy from NCBI for obitools pipeline into new TAXO folder 
    # note: it can't hurt to redownload as errors are possible when using older taxonomy databases
mkdir TAXO/
cd TAXO/
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
28. format NCBI taxonomy into obitools/ecoPCR format by running obitaxonomy.sh
obitaxonomy -t TAXO -d TAXO 
29. clean the ref database by running obigrep.ref1.sh (keeps only those ID'd to family level)
    # note: some of steps 29-32 may be unnecessary, but it at least serves to get it into obi format for ecotag
cd refs
obigrep -d ../TAXO/TAXO --require-rank=species \
  --require-rank=genus --require-rank=family ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln.C.fasta > ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.C.fasta
obigrep -d ../TAXO/TAXO --require-rank=species \
  --require-rank=genus --require-rank=family combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln.fasta > combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.fasta
30. dereplicate reference sequences by running obiuni.ref.sh (also adds full taxonomy to header)
obiuniq -d ../TAXO/TAXO \
  ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.C.fasta > ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.uni.C.fasta
31. ensure that the dereplicated sequences have a taxid at the family level by running obigrep.ref2.sh
obigrep -d ../TAXO2/TAXO2 --require-rank=family \
  ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.uni.C.fasta > ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.uni.cln.C.fasta
obigrep -d ../TAXO/TAXO --require-rank=family \
  combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.uni.fasta > combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.uni.cln.fasta
32. annotate reference database with unique IDs by running obiann.ref.sh
    # note: probably could run this one without a shell script
obiannotate --uniq-id ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.uni.cln.C.fasta > ncbi.bold.art.hum.mus.may23.pcr.pga.tax.uni.cln2.uni.cln.ann.C.fasta
obiannotate --uniq-id combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.uni.cln.fasta > combined.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.uni.cln.ann.fasta
obiannotate --uniq-id NJ_plus_exotics_local.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.uni.cln.fasta > NJ_plus_exotics_local.fish.reference.may23.pcr.sin.pga.tax.fix.uni.cln2.uni.cln.ann.fasta

