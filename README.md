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
for f in `ls -1 *.ali.cut.n21.fastq | sed 's/.ali.cut.n21.fastq//' `
do
seqtk seq -a ${f}.ali.cut.n21.fastq > ${f}.ali.cut.n21.fasta
done
```
# Remove primers, etc.
This is done using the obitools command ngsfilter. Use nano to make bash scripts to run the following code for each sequencing pool: f03_ngsfilter.X.sh, f03_ngsfilter.Y.sh, and f03_ngsfilter.Z.sh. Change the sequence pool code (3011 in the example below) in the file as needed. Run the jobs using sbatch called from within the rawdata folder. This script takes about a minute per file in its current form.
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

Merge files into one
```
cat *.ngs.fasta > merged.fasta
```
Dereplicate (takes about 30-40 min). Use nano to make a bash script (f04_mergeuni.sh) with the following code and run using sbatch.
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
obigrep -l 140 -L 190 -p 'count>=0' merged.uni.fasta > merged.uni.c0.140.190.fasta
obigrep -l 140 -L 190 -p 'count>=2' merged.uni.fasta > merged.uni.c2.140.190.fasta
obigrep -l 140 -L 190 -p 'count>=10' merged.uni.fasta > merged.uni.c10.140.190.fasta

```

```
obistat -c count merged.uni.c0.140.190.fasta |  sort -nk1 | head -20
```


18. ONLY DO THIS IF YOU DON't WANT TO USE SWARM TO CLUSTER. Run obiclean.sh to denoise (default settings? it is a simple form of clustering) - final file merged.uniq.c10.l140.l190.cln.fasta
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
optional: clustering sequences into MOTUs with SWARM (if doing this, don't run obiclean in # 18 above)

# Clustering into MOTUs

Move merged.uni.c0.140.190.fasta into a new folder called "cluster"
18b. shorten the uniq IDs (use a one-character prefix for it to work with the tallying script below)
```
obiannotate --seq-rank merged.uni.c10.140.190.fasta | obiannotate --set-identifier '"F_%d"%seq_rank' > merged.uni.c10.140.190.sht.fasta

obitab -o merged.uni.c10.140.190.sht.fasta > merged.uni.c10.140.190.sht.tab
```
Convert into vsearch format
```
obiannotate -k count merged.uni.c10.140.190.sht.fasta > merged.uni.c10.140.190.sht.vsc.fasta

sed -i -e 's/[[:space:]]count/;size/g' merged.uni.c10.140.190.sht.vsc.fasta
sed -i 's/;\+$//' merged.uni.c10.140.190.sht.vsc.fasta
```
23b. remove the chimeras
Make a script file (f05_chimera.sh) and run it using sbatch.
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

vsearch --uchime_denovo merged.uni.c10.140.190.sht.vsc.srt.fasta --sizeout --nonchimeras merged.uni.c10.140.190.sht.vsc.srt.chi.fasta --chimeras merged.uni.c10.140.190.sht.vsc.srt.chi.fasta --uchimeout uchimeout.txt
```
24b. make it into a single-line fasta file
```
awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' merged.uni.c10.140.190.sht.vsc.srt.chi.fasta > merged.uni.c10.140.190.sht.vsc.srt.chi.sin.fasta
```

25b. create the MOTUs (-d = max distance in nucleotides for 2 sequences to be same; setting at 1 based on Swarm github site advice)
swarm -d 1 -f -z -t 1 -o merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_output -s merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_stats -w merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1.fasta merged.uni.c10.140.190.sht.vsc.srt.chi.sin.fasta

26b. create R script program to recount SWARM reads by OTU and sample.

```
#!/FILEPATHTOTHESCRIPT/Rscript

## Script for recount abundances of a dataset clustered using Swarm to a tabulated csv file
## The script will read two arguments from the command line: the input file name (output file from swarm)
## and the abundances tab file from obitab.
## The output will be a ".csv" file with recalculated abundances per MOTU and per sample.
## By Owen S. Wangensteen - Project Metabarpark  2017
## Modified by Mike Allen, May 15, 2023, so that it would run with his data.

args<-commandArgs(T)

if (length(args)<2) {message("Please enter a cluster list output file obtained from SWARM and a tabulated counts file from obitab.")} else
{  
  fileswarm <- args[1] #Must be an output file from swarm 
  filetab <- args[2] 
  outfile <-paste(fileswarm,".counts.csv",sep="")   
  
  num_char_id <- 9 # changed from 14 by Mike so that it would work with his data
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
  clusters <- strsplit(swarm_db," ")
  for (i in 1:total_swarms) for (j in 1:length(clusters[[i]])) if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))==";"){
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-1)
  }
  cluster_reads  <- NULL
  for (i in 1:total_swarms) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
  swarm_db_reduced <- swarm_db[cluster_reads>=min_reads]
  
  clusters <- strsplit(swarm_db_reduced," ")
  total_swarms_reduced <- length(swarm_db_reduced)
  
  id <- NULL
  for (i in 1:total_swarms_reduced) for (j in 1:length(clusters[[i]])) {
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,num_char_id)  
    id[i] <- clusters[[i]][1]
  }

  ### inserted by Mike to remove ";" from cluster and sequence names
  for (i in 1:total_swarms_reduced) { for (j in 1:length(clusters[[i]])) {
    if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))==";"){
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-1)
    }
    if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))=="s"){
      clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-2)
    }
    if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))=="i"){
      clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-3)
    }
    if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))=="z"){
      clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-4)
    }
    if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))=="e"){
      clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-5)
    }
    if (substr(clusters[[i]][[j]],nchar(clusters[[i]][[j]]),nchar(clusters[[i]][[j]]))=="="){
      clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,nchar(clusters[[i]][[j]])-6)
    }
  }
    id[i] <- clusters[[i]][1]
  }
  ###

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
}
```

```
chmod +x f06_owi_recount_swarm

# change the /user/bin/Rscript at the top to your env URL (e.g., /home/USERNAME/miniconda3/envs/ren423/bin/Rscript)
# also need to change  num_char_id to 9 (from 14) and change the separator character for the output file to a "," 
# most importantly, a whole chunk of code needs to be inserted into the file to deal with cluster names (see Mike's version)
```
27b. recount the SWARM samples by MOTU (also excludes singletons - seqs with total read count = 1 across samples)
    # note: the .tab file was created in step 19b above.
```
f06_owi_recount_swarm merged.uni.c10.140.190.sht.vsc.srt.chi.sin.sw1_output merged.uni.c10.140.190.sht.tab
```
28b. fix header of the fasta file (add space between MOTU name and the size variable)
cp merged.uni.c10.l130.L185.sht.srt.nochi.1line.swarm1.Z.fasta merged.uni.c10.l130.L185.sht.srt.nochi.1line.swarm1.fix.Z.fasta
sed -i 's/;size=/ size=/g' merged.uni.c10.l130.L185.sht.srt.nochi.1line.swarm1.fix.Z.fasta
cp merged.uni.c10.l75.L125.sht.srt.nochi.1line.swarm1.C.fasta merged.uni.c10.l75.L125.sht.srt.nochi.1line.swarm1.fix.C.fasta
sed -i 's/;size=/ size=/g' merged.uni.c10.l75.L125.sht.srt.nochi.1line.swarm1.fix.C.fasta
cp merged.uni.c10.l140.L190.sht.srt.nochi.1line.swarm1.fasta merged.uni.c10.l140.L190.sht.srt.nochi.1line.swarm1.fix.fasta
sed -i 's/;size=/ size=/g' merged.uni.c10.l140.L190.sht.srt.nochi.1line.swarm1.fix.fasta
29. remove singletons from the SWARM fasta file (not necessary, I think, as we already cut seqs w/ <10 reads)
obigrep -p ‘size>1’ merged.uni.c10.l130.L185.sht.srt.nochi.1line.swarm7.fix.Z.fasta > merged.uni.c10.l130.L185.sht.srt.nochi.1line.swarm7.fix.no1.Z.fasta
obigrep -p ‘size>1’ merged.uni.c10.l75.L125.sht.srt.nochi.1line.swarm5.fix.C.fasta > merged.uni.c10.l75.L125.sht.srt.nochi.1line.swarm5.fix.no1.C.fasta

