# fish23
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
# install swarm (I think install this one on a different conda env with python 3+)
conda install -c bioconda swarm
# install ecoprimers
conda install -c bioconda ecoprimers

# Demultiplex data


# Get data onto cluster
```
wget -r -nH --cut-dirs=2 https://htseq.princeton.edu/tmp/SOMECODE/
```

# 
