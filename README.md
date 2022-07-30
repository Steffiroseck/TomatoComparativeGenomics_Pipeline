# tomato_comparative_genomics_pipeline
# Contents
* Introduction
* Prerequisites and Installation
* Usage

# Introduction
This is a pipeline which performs comparative genomics analysis of tomato and its wild relatives based on Gene Ontologies.
The main objective behind creating this pipeline was to find out the genetic changes between the domestic variety of tamato (_Solanum.lycopersicum_), and its wild relatives (_Solanum.chilense, Solanum.pennellii, Solanum.pimpinellifolium, Solanum.sitiens_, and _Solanum.lycopersicoides_) to better understand how the wild varities withstand abiotic stresses like salinity and drought, and why the domestic variety does not.

The initial file for this pipeline is a file which has the gene ids, their start and end coordinates, gene descriptions, and gene ontologies of _S. chilense_. The pipeline takes this initial file and extracts the Gene Ontolgy (GO) IDs and tries to find the IDs which are related to salinity or drought tolerances using an R package. Next, the genes and their sequences are extracted for the GO IDs having salinity or drought GO terms. These gene sequences are then compared with the whole genome of the domestic species _S.lycopersicum_ to find out sequence differences. 

# Prerequisites and Installation
1. R packages - GO.db, dplyr
2. samtools v1.9
3. NCBI Blast v2.13.0
4. MAFFT v7.490 
5. Snp-sites
6. bcftools v1.6 
7. seqkit
8. Provean 

Please follow the below instructions if the R packages and tools used in this pipeline are not already installed on your system.

# 1. GO.db package

    start R (version "4.2") and enter:

        install.packages("dplyr")
        if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install("GO.db")

# 2. samtools

Download the most current version from the Samtools website (http://sourceforge.net/projects/samtools/files/samtools/) 
    
Unzip the file  
    
In Bash,

		tar xvjf samtools-X.X.tar.bz2  
        
Go into the newly created directory and compile the code by typing make: 
    
		cd samtools-X.X     
        	make     
        
Modify your .bashrc file so that when you type "samtools" it calls the program: 
    
		export PATH=$PATH:/directory/samtools-X.X
        
# 3. NCBI-BLAST v2.13.0

This is available at the NCBI ftp site (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/)

# 4. MAFFT v7.490

Download the pre-compiled packages from the website (https://mafft.cbrc.jp/alignment/software/mafft-7.490-with-extensions-src.tgz)
    
Untar the package  
    
 		gunzip -cd mafft-x.x-src.tgz | tar xfv -      
     	cd mafft-x.x/core/   
        
Compile and install 
    
     	make       
     	make install  
        
Modify the .bashrc file 
    
     	export PATH=$PATH:/directory/mafft-x.x/bin/     
        
# 5. snp-sites. 

This tool has a git repository and can be directory cloned from the repository

Clone from the git respository. Type in bash;  

		git clone https://github.com/sanger-pathogens/snp-sites.git   
        
Install  
    
    	cd snp-sites/      
    	autoreconf -i -f      
    	./configure      
    	make      
    	make install  
        
Modify the .bashrc file  
    
    	export PATH=$PATH:/directory/snp-sites/src/  
        
# 6. bcftools

Download the most current version from the bcftools website (https://github.com/samtools/bcftools/releases/download/)
    
Unzip the file  
    
In Bash,
    
		tar xvjf bcftools-X.X.tar.bz2  
        
Go into the newly created directory and compile the code by typing make: 
    
    	cd bcftools-X.X     
    	make     
        
Modify your .bashrc file so that when you type "bcftools" it calls the program: 
    
    	export PATH=$PATH:/directory/bcftools-X.X

# 7. seqkit

# 8. Provean

PROVEAN requires the CD-HIT 3.1.2 (or more recent) software, ncbi-nr (non-redundant) protein database and NCBI_BLAST. So make sure these are installed in the system.
Visit http://weizhong-lab.ucsd.edu/cd-hit/download.php for installation guidelines for CD-HIT.
Current version of nr database is available at the NCBI ftp site ftp://ftp.ncbi.nih.gov/blast/db/

Download the Provean source code from http://sourceforge.net/projects/provean/

In bash,

		tar zxvf provean-1.1.5.tar.gz
		cd provean-1.1.5
		./configure --prefix=/path/to/where/provean/to/be/installed
		make
		make install
	

Alternatively, if your system has Anaconda or Miniconda installed, a conda installation of some of the softwares can be done. 

   		conda install -c bioconda samtools
   		conda install -c bioconda bcftools
   		conda install -c bioconda mafft
   		conda install -c bioconda seqkit
