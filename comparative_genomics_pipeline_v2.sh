#!/bash/bin

# path to the working directory
  wd=/home/steffi

# Paths to the softwares in the pipeline in the Rambo server
  blast=/opt/bix/blast+/2.13.0/bin
  seqkit=$wd/anaconda3/bin/seqkit
  mafft=$wd/installations/mafft-linux64/mafftdir/bin/mafft
  nrdb=$wd/nr_2.4/nr
  psiblast=/opt/bix/blast+/2.4.0/bin/psiblast
  cdhit=/opt/bix/cdhit/4.6.1/cd-hit
  blastdbcmd=/opt/bix/blast+/2.4.0/bin/blastdbcmd
  provean=/opt/bix/provean/1.1.5/src/provean


# Assuming that the initial file for analysis is in /home/steffi which is the present working directory.
# Create directories to store the results
  mkdir $wd/Results
  
# Create subdirectories to store results of different species comparisons
  mkdir Results/chilense_lycopersicum
  mkdir Results/chilense_lycopersicoides
  mkdir Results/chilense_pennellii
  mkdir Results/chilense_pimpinellifolium
  mkdir Results/chilense_sitiens

# Create variables to store the output directories
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

# Extract IDs and gene descriptions from the file
  #tr ";" "\n" < $wd/chilense_augustus_with_description.gff3 | grep -E '(ID|Description=)' | cut -d '=' -f 2 | paste - - |  tr ' ' '_' > $wd/chilense.id.description

# Extract the chilenseID , CDS, CDS start, CDS end and  Gene Ontology IDs and write to a file.
# The 9th column has the GO IDs, which are separated using ";". These GO IDs should be splitted.
# Split 9th column based on ;. This will make GO IDs in a separate column separated by tab and print only GO IDS
# Split the comma separated values in each column to multiple rows
# split the IDs based on "=" and ","
  cat $wd/chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1,$3,$4,$5, $9}' | sed 's/[; ]\+/\t/g' | awk '{for (i=1;i<=NF;i++){if ($i ~/Ontology_id/) {print $i}}}' | tr '=' '\n' | tr ',' '\n' > $wd/Results/Chilense.GO.IDs

# Now, remove .ontology_id. from the  file
  sed -i '/Ontology_id/d' $wd/Results/Chilense.GO.IDs

# Next, For all the GO IDs extracted from S.chilense, the GO Terms are extrcated. 
# R code to extract GO terms for the GO IDs extracted.
# This r script will save the GO IDs with salt/drought keywords to chilense.GO.terms.salt.drought file
  Rscript $wd/TomatoComparativeGenomics_Pipeline/extract_GO_terms.R

# Extract just the important GO IDs, so remove the GO terms from the R output file
  awk '{print $1}' $wd/Results/chilense.GO.terms.salt.drought | sed '1d' > $wd/Results/chilense.GO.terms.salt.drought.IDs
  
# Extract the gene IDs corresponding to each of the GO IDs related to salt/drought keywords
  grep -f $wd/Results/chilense.GO.terms.salt.drought.IDs <$wd/chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1}' > $wd/Results/chilense.GO.terms.salt.drought.IDs.genes

# Extract both the nucleotide and amino acid sequences. The sequences are extracted from the S.chilense genome coding sequence (augustus.with.hints.codingseq) and protein sequence (augustus.with.hints.filtered.aa.fasta) files which are stored in /home/steffi/ folder
# Nucleotide sequence extraction
  seqkit grep -f $wd/Results/chilense.GO.terms.salt.drought.IDs.genes $wd/augustus.with.hints.codingseq > $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta

# Save the sequences in a single line for tree generation
  grep -v ">" $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta|tr '\n' ' ' | sed -e 's/ //g' |sed 's/^/>S.chilense\n/'> $wd/Results/chilense.tree.fasta

# Amino acid sequence extraction
  seqkit grep -f $wd/Results/chilense.GO.terms.salt.drought.IDs.genes $wd/augustus.with.hints.filtered.aa.fasta > $wd/Results/chilense.GO.terms.salt.drought.IDs.aa.fa

# Merge S.chilense and other species proteins file
  cat augustus.with.hints.filtered.aa.fasta $lycopersicum_protein > chilense_lycopersicum_genome_aa.fa
  cat augustus.with.hints.filtered.aa.fasta $lycopersicoides_protein > chilense_lycopersicoides_genome_aa.fa
  cat augustus.with.hints.filtered.aa.fasta $pennellii_protein > chilense_pennellii_genome_aa.fa
  cat augustus.with.hints.filtered.aa.fasta $pimpinellifolium_protein > chilense_pimpinellifolium_genome_aa.fa
  cat augustus.with.hints.filtered.aa.fasta $sitiens_protein > chilense_sitiens_genome_aa.fa

1. Sequence comparison between S.chilense & S.lycopersicum

# Blastn is used for performing the sequence similarity search. Only the query coverage per HSP greater than 90% is reported.
# Make the Blast database for the Blastn analysis
  $blast/makeblastdb -in $lycopersicum -dbtype 'nucl'
  $blast/makeblastdb -in $pennellii -dbtype 'nucl'
  $blast/makeblastdb -in $lycopersicoides -dbtype 'nucl'
  $blast/makeblastdb -in $pimpinellifolium -dbtype 'nucl'
  $blast/makeblastdb -in $sitiens -dbtype 'nucl'

# Run Blastn
  $blast/blastn -db $lycopersicum -query $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cl/chilense.lycopersicum.blastn
  $blast/blastn -db $pennellii -query $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cp/chilense.pennellii.blastn
  $blast/blastn -db $lycopersicoides -query $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cly/chilense.lycopersicoides.blastn
  $blast/blastn -db $pimpinellifolium -query $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cpi/chilense.pimpinellifolium.blastn
  $blast/blastn -db $sitiens -query $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta -max_target_seqs 1 -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cs/chilense.sitiens.blastn


# Extract IDs of the wild and domestic species
  awk '{print $1}' $cl/chilense.lycopersicum.blastn > $cl/chilenseID_lycopersicum
  awk '{print $2}' $cl/chilense.lycopersicum.blastn > $cl/lycopersicumID
  awk '{print $1}' $cp/chilense.pennellii.blastn > $cp/chilenseID_pennellii
  awk '{print $2}' $cp/chilense.pennellii.blastn > $cp/pennelliiID
  awk '{print $1}' $cly/chilense.lycopersicoides.blastn > $cly/chilenseID_lycopersicoides
  awk '{print $2}' $cly/chilense.lycopersicoides.blastn > $cly/lycopersicoidesID
  awk '{print $1}' $cpi/chilense.pimpinellifolium.blastn > $cpi/chilenseID_pimpinellifolium
  awk '{print $2}' $cpi/chilense.pimpinellifolium.blastn > $cpi/pimpinellifoliumID
  awk '{print $1}' $cs/chilense.sitiens.blastn > $cs/chilenseID_sitiens
  awk '{print $2}' $cs/chilense.sitiens.blastn > $cs/sitiensID


# Get total number of unique ids from the Blastn results
  #sort $cl/chilenseID_lycopersicum | uniq | wc -l 
  #sort $cl/lycopersicumID | uniq | wc -l 

# Find chilenseIDs which had no mappings in S. lycopersicum after Blastn analysis
  #awk 'FNR==NR{a[$0]=1;next}!($0 in a)' $cl/chilense.lycopersicum.blastn Results/chilense.GO.terms.regions > $cl/chilenseIDs.not.present.blastn 

#get unique gene counts
  #sort $cl/chilenseIDs.not.present.blastn | uniq > $cl/chilenseIDs.not.present.blastn.unique

# Extract Query and Reference gene IDs from the Blastn results in each line to a separate file
  cd $cl/
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.lycopersicum.blastn
  cd $cp
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.pennellii.blastn
  cd $cly
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.lycopersicoides.blastn
  cd $cpi
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.pimpinellifolium.blastn
  cd $cs
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.sitiens.blastn


# Extract the aminoacid sequences (7 minutes to finish)
  for i in $cl/*.txt $cp/*.txt $cly/*.txt $cpi/*.txt $cs/*.txt; 
  do 
  seqkit grep -f ${i} $wd/chilense_lycopersicum_genome_aa.fa $wd/chilense_pennellii_genome_aa.fa $wd/chilense_lycopersicoides_genome_aa.fa $wd/chilense_pimpinellifolium_genome_aa.fa $wd/chilense_sitiens_genome_aa.fa > ${i%.txt*}.aminoacid.fa; 
  done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
  for i in $cl/*.aminoacid.fa $cp/*.aminoacid.fa $cly/*.aminoacid.fa $cpi/*.aminoacid.fa $cs/*.aminoacid.fa; 
  do 
  awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; 
  done

# Store the S.chilense gene sequence from the above file to a separate file. (This is used for variant effect prediction analysis)  
  for i in $cl/*.new.fa $cp/*.new.fa $cly/*.new.fa $cpi/*.new.fa $cs/*.new.fa; 
  do 
  awk '/^>/{if(N)exit;++N;} {print;}' ${i} > ${i%.new.fa*}.proveaninput.fa; 
  done
  
# Now perform MAFFT MSA (7 minutes to finish)
  for i in $cl/*.new.fa $cp/*.new.fa $cly/*.new.fa $cpi/*.new.fa $cs/*.new.fa; 
  do 
  mafft ${i} > ${i%.new.fa*}.mafft_out.fa; 
  done

# Extract the positions of variants, insertions, and deletions from the MSA file. Remove thw whitespaces from file so that output is compatible for running PROVEAN. 
  for i in $cl/*.mafft_out.fa $cp/*.mafft_out.fa $cly/*.mafft_out.fa $cpi/*.mafft_out.fa $cs/*.mafft_out.fa; 
  do 
  echo ${i} processing....
  python $wd/TomatoComparativeGenomics_Pipeline/extract_variants_indels_from_MSA_v2.py -f ${i} > ${i%.mafft_out.*}.var
  done

# Run PROVEAN analysis
  for i in $cl/*.proveaninput.fa $cp/*.proveaninput.fa $cly/*.proveaninput.fa $cpi/*.proveaninput.fa $cs/*.proveaninput.fa;
  do 
  cmd="$provean -q ${i} -d $nrdb -v ${i%.proveaninput.fa*}.var --psiblast $psiblast --cdhit $cdhit --blastdbcmd $blastdbcmd --num_threads 40 > ${i%.proveaninput.fa*}.proveanout"
  echo ${cmd}
  eval ${cmd}
  done
  cd

# Extract only the deleterious variants from the Provean output.
# Deleterious variants are the ones with PROVEAN score < -2.5.
# We also remove empty files which doesnt have any deleterious variants
  for i in $cl/*.proveanout $cp/*.proveanout $cly/*.proveanout $cpi/*.proveanout $cs/*.proveanout;
  do
  awk '/# VARIATION/{p=1}p' ${i} | awk -F'\t' '$2 < -2.5' > ${i%.proveanout*}.deleterious;
  done

# Provean analysis might generate empty files for the genes with no variants between the comparisons. You  can remove them by typing;
  find  -type f -empty -delete
  #find  -type f -name '*.deleterious' -empty | wc -l
  #find  -type f -name '*.var' -empty | wc -l
  
# Write all the deleterious variants with their number sorted to a file
  find $cl -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' > $cl/chilense_lycopersicum_deleteriousgenes
  find $cp -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' > $cp/chilense_pennellii_deleteriousgenes
  find $cly -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' > $cly/chilense_lycopersicoides_deleteriousgenes
  find $cpi -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' > $cpi/chilense_pimpinellifolium_deleteriousgenes
  find $cs -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' > $cs/chilense_sitiens_deleteriousgenes 
  
# Get the Top 5 genes with the highest number of deleterious variants
  find $cl -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' | head -5 > $cl/Top5_deleterious_genes
  find $cp -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' | head -5 > $cp/Top5_deleterious_genes
  find $cly -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' | head -5 > $cly/Top5_deleterious_genes
  find $cpi -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' | head -5 > $cpi/Top5_deleterious_genes
  find $cs -type f -name "*.deleterious" | xargs wc -l | sort -rn | grep -v ' total$' | head -5 > $cs/Top5_deleterious_genes  
  
# Finding the deleterious genes that are shared between salt and drought resistant and sensitive species
  for i in $cl/chilense_lycopersicum_deleteriousgenes $cp/chilense_pennellii_deleteriousgenes $cly/chilense_lycopersicoides_deleteriousgenes $cpi/chilense_pimpinellifolium_deleteriousgenes $cs/chilense_sitiens_deleteriousgenes;
  do
  sed 's|.*/\(.*\)_.*|\1|' $i > ${i%_deleteriousgenes*}_deleteriousgeneIDs
  done
  
  cat $cl/chilense_lycopersicum_deleteriousgeneIDs Results/chilense_lycopersicoides/chilense_lycopersicoides_deleteriousgeneIDs Results/chilense_pennellii/chilense_pennellii_deleteriousgeneIDs Results/chilense_pimpinellifolium/chilense_pimpinellifolium_deleteriousgeneIDs Results/chilense_sitiens/chilense_sitiens_deleteriousgeneIDs |sort |uniq -c |sed -n -e 's/^ *5 \(.*\)/\1/p' > $wd/Results/common.deleteriousgenes
  
  chilense_description=$wd/chilense.id.description
  grep -i -f $wd/Results/common.deleteriousgenes $chilense_description > $wd/Results/common.deleteriousgenes.descriptions
 
# Generate an Upset plot to get the gene IDs intersections between different species in the study. 
# A data table needs to be created with all the S.chilense IDs matched to the different species in the study. Once the data table is created, the R script to generate upset plot is saved as a separate script.
  paste $wd/Results/chilense_lycopersicum/chilenseID_lycopersicum $wd/Results/chilense_lycopersicoides/chilenseID_lycopersicoides $wd/Results/chilense_pennellii/chilenseID_pennellii $wd/Results/chilense_pimpinellifolium/chilenseID_pimpinellifolium $wd/Results/chilense_sitiens/chilenseID_sitiens > $wd/Results/GO.IDs.all.species

# Add column names
  sed  -i '1i chilense_lycopersicum\tchilense_lycopersicoides\tchilense_pennellii\tchilense_pimpinellifolium\tchilense_sitiens' $wd/Results/GO.IDs.all.species

# Run the Rscript to generate the upset plot
  Rscript $wd/TomatoComparativeGenomics_Pipeline/create_upsetplot.R

# Generate a phylogenetic tree of all the Blast hits in lycopersicum, lycopersicoides, pennellii, pimpinellifolium, and sitiens against chilense GO terms
# Extract the fasta sequences of the IDs in the Blast hits
  seqkit grep -f $cl/lycopersicumID $lycopersicum | awk 'BEGIN {print ">S.lycopersicum"};!/^>/{print}' > $cl/lycopersicumID.fasta
  seqkit grep -f $cly/lycopersicoidesID $lycopersicoides | awk 'BEGIN {print ">S.lycopersicoides"};!/^>/{print}' > $cly/lycopersicoidesID.fasta
  seqkit grep -f $cp/pennelliiID $pennellii | awk 'BEGIN {print ">S.pennelli"};!/^>/{print}' > $cp/pennelliiID.fasta
  seqkit grep -f $cpi/pimpinellifoliumID $pimpinellifolium | awk 'BEGIN {print ">S.pimpinellifolium"};!/^>/{print}' > $cpi/pimpinellifoliumID.fasta
  seqkit grep -f $cs/sitiensID $sitiens | awk 'BEGIN {print ">S.sitiens"};!/^>/{print}' > $cs/sitiensID.fasta
  
 # Concatenate the above fasta sequences into a single multifasta file.
 # This is for running the Multiple Sequence Alignment
   cat $wd/Results/chilense.tree.fasta $cl/lycopersicumID.fasta $cly/lycopersicoidesID.fasta $cp/pennelliiID.fasta $cpi/pimpinellifoliumID.fasta $cs/sitiensID.fasta > $wd/Results/LLPPS.fasta
 # count the number of sequences in the Multi fasta merged file. It should be same as the number of genomes in the study
   grep ">" $wd/Results/LLPPS.fasta
   
 # Run the Multiple Sequence Alignment program
   mafft $wd/Results/LLPPS.fasta > $wd/Results/LLPPS.mafft.out
   
 # Generate the Phylogenetic tree
   raxmlHPC -d -p 12345 -m GTRGAMMAI -s $wd/Results/LLPPS.mafft.out -n $wd/Results/LLPPS
   
 # Visualise the tree in any of the tree viewing softwares like Figtree or iTOL.
