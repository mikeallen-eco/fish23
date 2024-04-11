# fish23
fish metabarcoding bioinformatics workflow

# Set up software

1. Connected to the Amarel cluster and start an interactive session so you don't use the login nodes  
Connect to the VPN if you are off campus (AnyConnect software).

2. Go to the shell in OnDemand, the terminal if on a Mac, or MobaXterm (or Putty or similar), log into the cluster.

```
ssh YOURNETID@amarel.rutgers.edu
```
3. Start an interactive session.
```
srun -p main -N 1 -c 2 -n 1 -t 05:00:00 --pty /bin/bash  # type that into a terminal/console
```

3. Download miniconda if you don't have it installed already (for managing software versions)
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda -V # check that miniconda is installed

4. create a new conda environment that runs python v 2.7
conda create --name obi2 python=2.7
# activate the obi2 environment
conda activate obi2

4. install obitools, ecopcr, swarm, and ecoprimers
# install obitools version 2
conda install -c bioconda obitools
# install ecopcr (also requires python 2.7)
conda install -c bioconda ecopcr
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
