# tomato_comparative_genomics_pipeline
This is a pipeline which performs comparative genomics of tomato and its wild relatives based on Gene Ontologies.
The main objective behind creating this pipeline was to find out the genetic changes between domestic variety of tamato (S.lycopersicum), and its wild relatives (S.chilense, S.pennellii, S.pimpinellifolium, S.sitiens, and S.lycopersicoides) to better understand how the wild varities withstand stresses like salinity and drought, and why the domestic variety does not.

The initial file for this analysis is a file which has the S. chilense gene ids, their start and end coordinates , gene descriptions, gene ontologies etc..This pipeline takes this initial file and extracts the Gene Ontolgy (GO) IDs and tries to find the IDs which are related to salinity or drought tolerances. Using these GO Ids, the genes and their sequences are extracted and compared with the domestic species S.lycopersicum. We also compare four other wild relatives of S.lycopersicum which are; S.pennellii, S.lycopersicoides, S.pimpinellifolium, and S.sitiens. 

# Tools used in this pipeline
1. R packages - GO.db, dplyr
2. samtools v1.9
3. NCBI Blast v2.13.0
4. MAFFT v7.490 (2021/Oct/30)
5. Snp-sites
6. bcftools v1.6 
7. seqkit

# If the R packages are not installed, follow the below commands;
start R (version "4.2") and enter:

    install.packages("dplyr")
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("GO.db")

If the softwares used in the pipeline are not already installed, please find the  instructions to install them below

1. samtools
    a. Download the most current version from the Samtools website (http://sourceforge.net/projects/samtools/files/samtools/) 
    
    b. Unzip the file  
    
    c. In Bash,
    
        tar xvjf samtools-X.X.tar.bz2  
        
    d. Go into the newly created directory and compile the code by typing make: 
    
        cd samtools-1.1     
        make     
        
    e. Modify your .bashrc file so that when you type "samtools" it calls the program: 
    
        export PATH=$PATH:/directory/samtools-X.X.XX 
        
2. NCBI-BLAST

4. MAFFT v7.490

    a. Download the pre-compiled packages from the website (https://mafft.cbrc.jp/alignment/software/mafft-7.490-with-extensions-src.tgz)
    
    b. Untar the package  
    
        gunzip -cd mafft-x.x-src.tgz | tar xfv -      
        cd mafft-x.x/core/   
        
    c. Compile and install 
    
        make       
        make install  
        
    d. Modify the .bashrc file 
    
        export PATH=$PATH:/directory/mafft-x.x/bin/     
        
4. snp-sites. This tool has a git repository and can be directory cloned from the repository

    a. Clone from the git respository. Type in bash;  
    
        git clone https://github.com/sanger-pathogens/snp-sites.git   
        
    b. Install  
    
        cd snp-sites/      
        autoreconf -i -f      
        ./configure      
        make      
        make install  
        
    c. Modify the .bashrc file  
    
        export PATH=$PATH:/directory/snp-sites/src/  
        
5. bcftools

6. seqkit

Alternatively, if your system has Anaconda or Miniconda installed, a conda installation of some of the softwares can be done. 

    conda install -c bioconda samtools
    conda install -c bioconda bcftools
    conda install -c bioconda mafft
    conda install -c bioconda seqkit
