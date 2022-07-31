The below script can be used to extract the sequence differences like inversions, gaps, breaks, duplications etc.. between the species comparisons. 
It uses a publicly available command-line software called MUMmer version 4.0.0. The gene IDs from S.lycopersicum, S.lycopersicoides, S.pennellii, S.pimpinellifolium, and S.sitiens that showed matches in the Blastn analysis are extracted and fasta sequences are retrieved.
These sequences are then run against the S. chilense salt/drought tolerant gene sequences and differences are tabulated.

# Path to current working directory
wd=/home/steffi
# Path to MUMmer software is stored in a variable
mummer_path="/opt/bix/mummer/4.0.0/bin/nucmer"

# Make a new directory to store the results from MUMmer analysis
mkdir $wd/Results/MUMmer_GOanalysis

# Extract fasta sequences
seqkit grep -f Results/chilense_lycopersicum/lycopersicumID /home/steffi/genomes/tomato/ITAG4.1_CDS.fasta > Results/nucmer_GOanalysis/lycopersicum.GOterms.matches.fasta
 
# Run nucmer alignment
/opt/bix/mummer/4.0.0/bin/nucmer --mum --minmatch=50 --mincluster=100 --prefix=Results/nucmer_GOanalysis/chilense.lycopersicum.nucmer --threads=8 Results/chilense.GO.terms.salt.drought.IDs.fasta Results/nucmer_GOanalysis/lycopersicum.GOterms.matches.fasta
 
# Filter the alignments
/opt/bix/mummer/4.0.0/bin/delta-filter -1 Results/nucmer_GOanalysis/chilense.lycopersicum.nucmer.delta > Results/nucmer_GOanalysis/chilense.lycopersicum.nucmer.delta.filter
 
# Run show-diff
/opt/bix/mummer/4.0.0/bin/show-diff -H -q Results/nucmer_GOanalysis/chilense.lycopersicum.nucmer.delta > Results/nucmer_GOanalysis/chilense.lycopersicum.nucmer.delta.diff #No inversions
