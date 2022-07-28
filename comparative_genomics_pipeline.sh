# path to the working directory
  wd=/home/steffi

# Paths to the softwares in the pipeline in the Rambo server
  samtools=/opt/bix/samtools/1.9/samtools
  blast=$wd/installations/ncbi-blast-2.13.0+/bin/blast
  seqkit=$wd/anaconda3/bin/seqkit
  bcftools=/opt/bix/bcftools/1.6/bcftools
  mafft=$wd/installations/mafft-linux64/mafftdir/bin/mafft
  snp-sites=$wd/installations/snp-sites/src/snp-sites


# Assuming that the initial file for analysis is in /home/steffi which is the present working directory.
# Create directories to store the results
mkdir $wd/Results
# Create subdirectories to store results of different species comparisons
mkdir Results/chilense_lycopersicum
mkdir Results/chilense_lycopersicoides
mkdir Results/chilense_pennellii
mkdir Results/chilense_pimpinellifolium
mkdir Results/chilense_sitiens

cl=Results/chilense_lycopersicum
cly=Results/chilense_lycopersicoides
cp=Results/chilense_pennellii
cpi=Results/chilense_pimpinellifolium
cs=Results/chilense_sitiens

# Store the genome CDS sequences and protein sequences of different species under the study to a variable
lycopersicum=$wd/genomes/tomato/ITAG4.1_CDS.fasta
lycopersicum_protein=$wd/genomes/tomato/ITAG4.1_proteins.fasta
lycopersicoides=$wd/genomes/Slycopersicoides/SlydLA2951_v2.0_cds.fasta
lycopersicoides_protein=$wd/genomes/Slycopersicoides/SlydLA2951_v2.0_protein.fasta
pennellii=$wd/genomes/Spennellii/Spenn-v2-cds-annot.fa
pennellii_protein=$wd/genomes/Spennellii/Spenn-v2-aa-annot.fa
pimpinellifolium=$wd/genomes/Spimpinellifolium/Spimpinellifolium_genome.CDS.fa
pimpinellifolium_protein=$wd/genomes/Spimpinellifolium/Spimpinellifolium_genome.protein.fa
sitiens=$wd/genomes/Ssitiens/augustus.hints.mrna
sitiens_protein=$wd/genomes/Ssitiens/sitiens_augustus.hints.aa


# Extract the chilenseID , CDS, CDS start, CDS end and  Gene Ontology IDs and write to a file.
# The 9th column has the GO IDs, which are separated using ";". These GO IDs should be splitted.
cat chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$3,$4,$5, $9}' > Results/Chilense.augustus.interpro.GO.IDs  

# Split 9th column based on ;. This will make GO IDs in a separate column separated by tab
sed 's/[; ]\+/\t/g' Results/Chilense.augustus.interpro.GO.IDs > Results/Chilense.augustus.interpro.GO.IDs.split

# Print only GO IDS and write it to a file. The output file will now look like each row with GO IDs corresponding to it separated by comma.
awk '{for (i=1;i<=NF;i++){if ($i ~/Ontology_id/) {print $i}}}' Results/Chilense.augustus.interpro.GO.IDs.split > Results/Chilense.GO.IDs

# Split the comma separated values in each column to multiple rows
# split the IDs based on "="
tr '=' '\n' < Results/Chilense.GO.IDs > Results/Chilense.GO.IDs.1 
# split based on .,.
tr ',' '\n' < Results/Chilense.GO.IDs.1 > Results/Chilense.GO.IDs.final

# Now, remove .ontology_id. from the  file
sed -i '/Ontology_id/d' Results/Chilense.GO.IDs.final #89582 GO IDs

# Next, For all the GO IDs extracted from S.chilense, the GO Terms are extrcated. This is done using GO.db (https://doi.org/doi:10.18129/B9.bioc.GO.db) package in R. 
# R code to extract GO terms for the GO IDs extracted
Rscript scripts/extract_GO_terms.R


# Extract ID, Gene start and stop from the Interpro.gff3 file for the GO IDs in the chilense.GO.terms.salt.drought
awk '{print $1}' Results/chilense.GO.terms.salt.drought | sed '1d' > Results/chilense.GO.terms.salt.drought.IDs

grep "GO:1901002" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:1901002"}' > Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0009414" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0009414"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0009651" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0009651"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0071472" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0071472"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:1901000" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:1901000"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:1901001" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:1901001"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0016717" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0016717"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0042631" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0042631"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0042538" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0042538"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0009819" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0009819"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0006833" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0006833"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0015250" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0015250"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:1902584" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:1902584"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0080148" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0080148"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:2000070" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:2000070"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0009415" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0009415"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0050521" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0050521"}' >> Results/chilense.GO.terms.salt.drought.IDs.location
grep "GO:0050891" chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$4,$5, "GO:0050891"}' >> Results/chilense.GO.terms.salt.drought.IDs.location

# Extract sequences from file. Inorder to do this, the file with GO IDs should be in a correct format. To get it, run the code below;
#cat Results/chilense.GO.terms.salt.drought.IDs.location | awk 'BEGIN {OFS = "\t"};  { printf("%s:%d-%d\n",$1,$2,$3)}' > Results/regions.txt

# Extract full length gene sequences of S.chilense based on the IDs found to be involved in salt or drought tolerances
cat Results/chilense.GO.terms.salt.drought.IDs.location | awk 'BEGIN {OFS = "\t"};  { printf("%s\n",$1,$2,$3)}' > Results/chilense.GO.terms.regions

# Use samtools to extract both the nucleotide and amino acid sequences. The sequences are extracted from the S.chilense genome coding sequence (augustus.with.hints.codingseq) and protein sequence (augustus.with.hints.filtered.aa.fasta) files which are stored in /home/steffi/ folder
# Nucleotide sequence extraction
#/opt/bix/samtools/1.9/samtools faidx augustus.with.hints.codingseq -r Results/chilense.GOtermsID.regions > Results/chilense.GO.terms.salt.drought.IDs.fasta
seqkit grep -f Results/chilense.GO.terms.regions augustus.with.hints.codingseq > Results/chilense.GO.terms.salt.drought.IDs.fasta

# Amino acid sequence extraction
#/opt/bix/samtools/1.9/samtools faidx augustus.with.hints.filtered.aa.fasta -r Results/chilense.GO.terms.regions > Results/chilense.GO.terms.salt.drought.IDs.aa.fa
seqkit grep -f Results/chilense.GO.terms.regions augustus.with.hints.filtered.aa.fasta > Results/chilense.GO.terms.salt.drought.IDs.aa.fa

# Merge S.chilense and other species proteins file
cat augustus.with.hints.filtered.aa.fasta /home/steffi/genomes/tomato/ITAG4.1_proteins.fasta > chilense_lycopersicum_genome_aa.fa
cat augustus.with.hints.filtered.aa.fasta $lycopersicoides_protein > chilense_lycopersicoides_genome_aa.fa
cat augustus.with.hints.filtered.aa.fasta $pennellii_protein > chilense_pennellii_genome_aa.fa
cat augustus.with.hints.filtered.aa.fasta $pimpinellifolium_protein > chilense_pimpinellifolium_genome_aa.fa
cat augustus.with.hints.filtered.aa.fasta $sitiens_protein > chilense_sitiens_genome_aa.fa

######################################################################################################################################################### 1. Sequence comparison between S.chilense & S.lycopersicum
# Blastn is used for performing the sequence similarity search. Only the query coverage per HSP greater than 90% is reported.

# Make the Blast database for the Blastn analysis
makeblastdb -in $lycopersicum -dbtype 'nucl'

# Run Blastn
blastn -db /home/steffi/genomes/tomato/ITAG4.1_CDS.fasta -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cl/chilense.lycopersicum.blastn #198 matches

# Extract just the chilenseIDs and lycopersicum IDs to a file
awk '{print $1}' $cl/chilense.lycopersicum.blastn > $cl/chilenseID_lycopersicum
awk '{print $2}' $cl/chilense.lycopersicum.blastn > $cl/lycopersicumID

# Get total number of unique ids from the Blastn result
sort $cl/chilenseIDs | uniq | wc -l #112
sort $cl/lycopersicumIDs | uniq | wc -l #111

# Find chilenseIDs which had no mappings in S. lycopersicum after Blastn analysis
awk 'FNR==NR{a[$0]=1;next}!($0 in a)' $cl/chilense.lycopersicum.blastn Results/chilense.GO.terms.regions > $cl/chilenseIDs.not.present.blastn #68

#get unique gene counts
sort $cl/chilenseIDs.not.present.blastn | uniq > $cl/chilenseIDs.not.present.blastn.unique #55
# Get gene descriptions for the unique genes in chilense

# Extract S.chilense and S.lycopersicum gene IDs from the Blastn results in each line to a separate file
cd $cl/
awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.lycopersicum.blastn

# Extract the aminoacid sequences
#for i in *.txt; do /opt/bix/samtools/1.9/samtools faidx chilense_lycopersicum_genome_aa.fa -r ${i} > ${i%.txt*}.aminoacid.fa; done
for i in *.txt; do seqkit grep -f ${i} /home/steffi/chilense_lycopersicum_genome_aa.fa > ${i%.txt*}.aminoacid.fa; done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
for i in *.aminoacid.fa ; do awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; done

# Now perform MAFFT MSA
for i in *.new.fa; do mafft ${i} > ${i%.new.fa*}.mafft_out.fa; done

# Run snp-sites on each mafft output file
for i in *.mafft_out.fa; do snp-sites -v -o ${i%.mafft_out.fa*}.snp-sites ${i} | echo ${i} completed; done #21 sequences had no aminoacid changes between chilense and lycopersicum

# Run the next step to get minoacid changes for PROVEAN analysis
for i in *.snp-sites; do /opt/bix/bcftools/1.6/bcftools query -f '%REF%POS%ALT\n' ${i} > ${i%.snp-sites*}.proveaninput; done
cd
# Run PROVEAN analysis
# Copy the proveaninput file contents for each file in the PROVEAN website's variants section. Copy the chilense amino acid gene sequence to the protein sequence section. Submit job 

###################################################################### RUN NUCMER BETWEEN GO TERMS in Chilense and lycopersicum
# copy the solIDS from blastn result againts chilense GO terms
#awk -F '\t' '{print $2}' chilense_GOtermsID_lycopersicum.blastn >solID_GOterms_matches

# Extract fasta sequences
# seqkit grep -f solID_GOterms_matches /home/steffi/genomes/tomato/ITAG4.1_CDS.fasta > solID_GOterms_matches.fasta
 
# Run nucmer alignment
#/opt/bix/mummer/4.0.0/bin/nucmer --mum --minmatch=50 --mincluster=100 --prefix=chilense_lycopersicum_nucmer --threads=8 chilense_GOtermsID.fasta solID_GOterms_matches.fasta
 
# Filter the alignments
#/opt/bix/mummer/4.0.0/bin/delta-filter -1 chilense_lycopersicum_nucmer.delta > chilense_lycopersicum_nucmer.delta.filter
 
# Run show-diff
#/opt/bix/mummer/4.0.0/bin/show-diff -H -q chilense_lycopersicum_nucmer.delta > chilense_lycopersicum_nucmer.delta.diff #No inversions

#######################################################################################################################################################

# 2. Sequence comparison between S.chilense and S.pennellii
# The GO terms fasta file of S.chilense is the query file for the analysis

# Blast similarity search of GO terms in S.chilense against S.pennellii
makeblastdb -in $pennellii -dbtype 'nucl'
blastn -db $pennellii -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cp/chilense.pennellii.blastn #265 matches

# Extract chilenseIDs and pennellii IDs to a file
awk '{print $1}' $cp/chilense.pennellii.blastn > $cp/chilenseID_pennellii
awk '{print $2}' $cp/chilense.pennellii.blastn > $cp/pennelliiID

# Extract chilense and sol IDS in each line to a separate file
cd $cp
awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.pennellii.blastn

# Extract the aminoacid sequences
for i in *.txt; do seqkit grep -f ${i} /home/steffi/chilense_pennellii_genome_aa.fa > ${i%.txt*}.aminoacid.fa; done
#for i in *.txt; do /opt/bix/samtools/1.9/samtools faidx /home/steffi/chilense_pennellii_genome_aa.fa -r ${i} > ${i%.txt*}.aminoacid.fa; done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
for i in *.aminoacid.fa ; do awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; done

# Now perform MAFFT MSA
for i in *.new.fa; do mafft ${i} > ${i%.new.fa*}.mafft_out.fa; done

# Run snp-sites on each mafft output file
for i in *.mafft_out.fa; do snp-sites -v -o ${i%.mafft_out.fa*}.snp-sites ${i} | echo ${i} completed; done #21 sequences had no aminoacid changes between chilense and lycopersicum

# Run the next step to get minoacid changes for PROVEAN analysis
for i in *.snp-sites; do /opt/bix/bcftools/1.6/bcftools query -f '%REF%POS%ALT\n' ${i} > ${i%.snp-sites*}.proveaninput; done
cd

#######################################################################################################################################################
# 3. Sequence comparison between S.chilense and S.lycopersicoides
# The GO terms fasta file of S.chilense is the query file for the analysis

# Blast similarity search of GO terms in S.chilense against S.pennellii
makeblastdb -in $lycopersicoides -dbtype 'nucl'
blastn -db $lycopersicoides -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cly/chilense.lycopersicoides.blastn #233 matches

# Extract chilenseIDs and pennellii IDs to a file
awk '{print $1}' $cly/chilense.lycopersicoides.blastn > $cly/chilenseID_lycopersicoides
awk '{print $2}' $cly/chilense.lycopersicoides.blastn > $cly/lycopersicoidesID

# Extract chilense and sol IDS in each line to a separate file
cd $cly
awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.lycopersicoides.blastn

# Extract the aminoacid sequences
for i in *.txt; do seqkit grep -f ${i} $wd/chilense_lycopersicoides_genome_aa.fa > ${i%.txt*}.aminoacid.fa; done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
for i in *.aminoacid.fa ; do awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; done

# Now perform MAFFT MSA
for i in *.new.fa; do mafft ${i} > ${i%.new.fa*}.mafft_out.fa; done

# Run snp-sites on each mafft output file
for i in *.mafft_out.fa; do snp-sites -v -o ${i%.mafft_out.fa*}.snp-sites ${i} | echo ${i} completed; done #21 sequences had no aminoacid changes between chilense and lycopersicum

# Run the next step to get minoacid changes for PROVEAN analysis
for i in *.snp-sites; do /opt/bix/bcftools/1.6/bcftools query -f '%REF%POS%ALT\n' ${i} > ${i%.snp-sites*}.proveaninput; done
cd

#######################################################################################################################################################
# 4. Sequence comparison between S.chilense and S.pimpinelifolium (LA1589)
# The GO terms fasta file of S.chilense is the query file for the analysis

# Blast similarity search of GO terms in S.chilense against S.pennellii
makeblastdb -in $pimpinellifolium -dbtype 'nucl'
blastn -db $pimpinellifolium -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cpi/chilense.pimpinellifolium.blastn #204 matches

# Extract chilenseIDs and pennellii IDs to a file
awk '{print $1}' $cpi/chilense.pimpinellifolium.blastn > $cpi/chilenseID_pimpinellifolium
awk '{print $2}' $cpi/chilense.pimpinellifolium.blastn > $cpi/pimpinellifoliumID

# Extract chilense and sol IDS in each line to a separate file
cd $cpi
awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.pimpinellifolium.blastn

# Extract the aminoacid sequences
for i in *.txt; do seqkit grep -f ${i} $wd/chilense_pimpinellifolium_genome_aa.fa > ${i%.txt*}.aminoacid.fa; done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
for i in *.aminoacid.fa ; do awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; done

# Now perform MAFFT MSA
for i in *.new.fa; do mafft ${i} > ${i%.new.fa*}.mafft_out.fa; done

# Run snp-sites on each mafft output file
for i in *.mafft_out.fa; do snp-sites -v -o ${i%.mafft_out.fa*}.snp-sites ${i} | echo ${i} completed; done #21 sequences had no aminoacid changes between chilense and lycopersicum

# Run the next step to get minoacid changes for PROVEAN analysis
for i in *.snp-sites; do /opt/bix/bcftools/1.6/bcftools query -f '%REF%POS%ALT\n' ${i} > ${i%.snp-sites*}.proveaninput; done
cd

######################################################################################################################################################## 4. Sequence comparison between S.chilense and S.sitiens
# The GO terms fasta file of S.chilense is the query file for the analysis

# Blast similarity search of GO terms in S.chilense against S.pennellii
makeblastdb -in $sitiens -dbtype 'nucl'
blastn -db $sitiens -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cs/chilense.sitiens.blastn #279 matches

# Extract chilenseIDs and pennellii IDs to a file
awk '{print $1}' $cs/chilense.sitiens.blastn > $cs/chilenseID_sitiens
awk '{print $2}' $cs/chilense.sitiens.blastn > $cs/sitiensID

# Extract chilense and sol IDS in each line to a separate file
cd $cs
awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.sitiens.blastn

# Extract the aminoacid sequences
for i in *.txt; do seqkit grep -f ${i} $wd/chilense_sitiens_genome_aa.fa > ${i%.txt*}.aminoacid.fa; done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
for i in *.aminoacid.fa ; do awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; done

# Now perform MAFFT MSA
for i in *.new.fa; do mafft ${i} > ${i%.new.fa*}.mafft_out.fa; done

# Run snp-sites on each mafft output file
for i in *.mafft_out.fa; do snp-sites -v -o ${i%.mafft_out.fa*}.snp-sites ${i} | echo ${i} completed; done #21 sequences had no aminoacid changes between chilense and lycopersicum

# Run the next step to get minoacid changes for PROVEAN analysis
for i in *.snp-sites; do /opt/bix/bcftools/1.6/bcftools query -f '%REF%POS%ALT\n' ${i} > ${i%.snp-sites*}.proveaninput; done
cd
