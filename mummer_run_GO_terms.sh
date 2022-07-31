The below script can be used to extract the sequence differences like inversions, gaps, breaks, duplications etc.. between the species comparisons. 
It uses a publicly available command-line software called MUMmer version 4.0.0. The gene IDs from S.lycopersicum, S.lycopersicoides, S.pennellii, S.pimpinellifolium, and S.sitiens that showed matches in the Blastn analysis are extracted and fasta sequences are retrieved.
These sequences are then run against the S. chilense salt/drought tolerant gene sequences and differences are tabulated.

# Path to current working directory
wd=/home/steffi

# Path to MUMmer software is stored in a variable
mummer_path="/opt/bix/mummer/4.0.0/bin"

# Path to whole genome files
lycopersicum=$wd/genomes/tomato/ITAG4.1_CDS.fasta
lycopersicoides=$wd/genomes/Slycopersicoides/SlydLA2951_v2.0_cds.fasta
pennellii=$wd/genomes/Spennellii/Spenn-v2-cds-annot.fa
pimpinellifolium=$wd/genomes/Spimpinellifolium/Spimpinellifolium_genome.CDS.fa
sitiens=$wd/genomes/Ssitiens/augustus.hints.mrna

# Make a new directory to store the results from MUMmer analysis and save the directory to a variable
mkdir $wd/Results/MUMmer_GOanalysis
output_dir=$wd/Results/MUMmer_GOanalysis

# Extract fasta sequences for the all the gene IDs from Blastn results.
seqkit grep -f $wd/Results/chilense_lycopersicum/lycopersicumID $lycopersicum > $output_dir/lycopersicum.GOterms.matches.fasta
seqkit grep -f $wd/Results/chilense_pennellii/pennelliiID $pennellii > $output_dir/pennellii.GOterms.matches.fasta
seqkit grep -f $wd/Results/chilense_lycopersicoides/lycopersicoidesID $lycopersicoides > $output_dir/lycopersicoides.GOterms.matches.fasta
seqkit grep -f $wd/Results/chilense_pimpinellifolium/pimpinellifoliumID $pimpinellifolium > $output_dir/pimpinellifolium.GOterms.matches.fasta
seqkit grep -f $wd/Results/chilense_sitiens/sitiensID $sitiens > $output_dir/sitiens.GOterms.matches.fasta

# Run nucmer alignment for each of the fasta file against S.chilense salt/drought GO terms fasta file
cd $output_dir/
for i in *.fasta; 
do $mummer_path/nucmer  --mum --minmatch=50 --mincluster=100 --prefix=${i%.GOterms.matches.fasta*}.chilense --threads=8 $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta ${i};
done

# Filter the alignments
for i in *.delta;
do $mummer_path/delta-filter -1 ${i} > ${i}.filter;
done
 
# Run show-diff to get the differences between the genome sequences
for i in *.delta;
do $mummer_path/show-diff -H -q ${i} > ${i}.diff;
done
