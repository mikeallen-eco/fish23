#!/FILEPATHTOTHESCRIPT/bin/Rscript

## Script for recount abundances of a dataset clustered using Swarm to a tabulated csv file
## The script will read two arguments from the command line: the input file name (output file from swarm)
## and the abundances tab file from obitab.
## The output will be a ".csv" file with recalculated abundances per MOTU and per sample.
## By Owen S. Wangensteen - Project Metabarpark  2017
## Modified by Mike Allen, so that it would run with his data.

args<-commandArgs(T)

if (length(args)<2) {message("Please enter a cluster list output file obtained from SWARM and a tabulated counts file from obitab.")} else
{  
  
  fileswarm <- args[1] #Must be an output file from swarm
  filetab <- args[2]
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
}