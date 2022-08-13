#!/bash/bin

# path to the working directory
  wd=/home/steffi

# Paths to the softwares in the pipeline in the Rambo server
  blast=/opt/bix/blast+/2.13.0/bin
  seqkit=$wd/anaconda3/bin/seqkit
  mafft=$wd/installations/mafft-linux64/mafftdir/bin/mafft
  nrdb=$wd/nr_2.4
  psiblast=/opt/bix/blast+/2.4.0/bin/psiblast
  cdhit=/opt/bix/cdhit/4.6.1/cd-hit
  blastdbcmd=/opt/bix/blast+/2.4.0/bin/blastdbcmd
  provean =/opt/bix/provean/1.1.5/src/provean


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
  Rscript $wd/scripts/extract_GO_terms.R

# Extract just the important GO IDs, so remove the GO terms from the R output file
  awk '{print $1}' $wd/Results/chilense.GO.terms.salt.drought | sed '1d' > $wd/Results/chilense.GO.terms.salt.drought.IDs
  
# Extract the gene IDs corresponding to each of the GO IDs related to salt/drought keywords
  grep -f $wd/Results/chilense.GO.terms.salt.drought.IDs <$wd/chilense_augustus_with_Interpro_filtered.gff3 | awk 'BEGIN {OFS = "\t"}; {print $1}' > $wd/Results/chilense.GO.terms.salt.drought.IDs.genes

# Extract both the nucleotide and amino acid sequences. The sequences are extracted from the S.chilense genome coding sequence (augustus.with.hints.codingseq) and protein sequence (augustus.with.hints.filtered.aa.fasta) files which are stored in /home/steffi/ folder
# Nucleotide sequence extraction
  seqkit grep -f $wd/Results/chilense.GO.terms.salt.drought.IDs.genes $wd/augustus.with.hints.codingseq > $wd/Results/chilense.GO.terms.salt.drought.IDs.fasta
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

# Run Blastn
  $blast/blastn -db /home/steffi/genomes/tomato/ITAG4.1_CDS.fasta -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cl/chilense.lycopersicum.blastn

# Extract just the chilenseIDs and lycopersicum IDs to a file
  awk '{print $1}' $cl/chilense.lycopersicum.blastn > $cl/chilenseID_lycopersicum
  awk '{print $2}' $cl/chilense.lycopersicum.blastn > $cl/lycopersicumID

# Get total number of unique ids from the Blastn result
  sort $cl/chilenseID_lycopersicum | uniq | wc -l 
  sort $cl/lycopersicumID | uniq | wc -l 

# Find chilenseIDs which had no mappings in S. lycopersicum after Blastn analysis
  awk 'FNR==NR{a[$0]=1;next}!($0 in a)' $cl/chilense.lycopersicum.blastn Results/chilense.GO.terms.regions > $cl/chilenseIDs.not.present.blastn 

#get unique gene counts
  sort $cl/chilenseIDs.not.present.blastn | uniq > $cl/chilenseIDs.not.present.blastn.unique

# Extract S.chilense and S.lycopersicum gene IDs from the Blastn results in each line to a separate file
  cd $cl/
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.lycopersicum.blastn

# Extract the aminoacid sequences
#for i in *.txt; do /opt/bix/samtools/1.9/samtools faidx chilense_lycopersicum_genome_aa.fa -r ${i} > ${i%.txt*}.aminoacid.fa; done
  for i in *.txt; 
  do 
  seqkit grep -f ${i} /home/steffi/chilense_lycopersicum_genome_aa.fa > ${i%.txt*}.aminoacid.fa; 
  done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
  for i in *.aminoacid.fa ; 
  do 
  awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; 
  done

# Store the S.chilense gene sequence from the above file to a separate file. (This is used for variant effect prediction analysis)  
  for i in *.new.fa; 
  do 
  awk '/^>/{if(N)exit;++N;} {print;}' ${i} > ${i%.new.fa*}.proveaninput.fa; 
  done
  
# Now perform MAFFT MSA
  for i in *.new.fa; 
  do 
  mafft ${i} > ${i%.new.fa*}.mafft_out.fa; 
  done

# Extract the positions of variants, insertions, and deletions from the MSA file. Remove thw whitespaces from file so that output is compatible for running PROVEAN. 
  for i in *.mafft_out.fa; 
  do 
  echo ${i} processing....
  python $wd/scripts/extract_variants_indels_from_MSA_v2.py ${i} > ${i%.mafft_out.*}.var
  done

# Run PROVEAN analysis
  for i in *.proveaninput.fa;
  do 
  cmd="$provean -q ${i} -d $nrdb -v ${i%.proveaninput.fa*}.var --psiblast $psiblast --cdhit $cdhit --blastdbcmd $blastdbcmd --num_threads 40 > ${i%.proveaninput.fa*}.proveanout"
  echo ${cmd}
  eval ${cmd}
  done
  cd

2. Sequence comparison between S.chilense and S.pennellii
# The GO terms fasta file of S.chilense is the query file for the analysis
# Blast similarity search of GO terms in S.chilense against S.pennellii
  $blast/makeblastdb -in $pennellii -dbtype 'nucl'
  $blast/blastn -db $pennellii -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cp/chilense.pennellii.blastn

# Extract chilenseIDs and pennellii IDs to a file
  awk '{print $1}' $cp/chilense.pennellii.blastn > $cp/chilenseID_pennellii
  awk '{print $2}' $cp/chilense.pennellii.blastn > $cp/pennelliiID

# Extract chilense and pennellii IDS in each line to a separate file
  cd $cp
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.pennellii.blastn

# Extract the aminoacid sequences
  for i in *.txt; 
  do 
  seqkit grep -f ${i} /home/steffi/chilense_pennellii_genome_aa.fa > ${i%.txt*}.aminoacid.fa; 
  done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
  for i in *.aminoacid.fa ; 
  do 
  awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; 
  done

# Store the S.chilense gene sequence from the above file to a separate file. (This is used for variant effect prediction analysis)  
  for i in *.new.fa; 
  do 
  awk '/^>/{if(N)exit;++N;} {print;}' ${i} > ${i%.new.fa*}.proveaninput.fa; 
  done
  
# Now perform MAFFT MSA
  for i in *.new.fa; 
  do 
  mafft ${i} > ${i%.new.fa*}.mafft_out.fa; 
  done

# Extract the positions of variants, insertions, and deletions from the MSA file. Remove thw whitespaces from file so that output is compatible for running PROVEAN. 
  for i in *.mafft_out.fa; 
  do 
  echo ${i} processing....
  python $wd/scripts/extract_variants_indels_from_MSA_v2.py ${i} > ${i%.mafft_out.*}.var
  done

# Run PROVEAN analysis
  for i in *.proveaninput.fa;
  do 
  cmd="$provean -q ${i} -d $nrdb -v ${i%.proveaninput.fa*}.var --psiblast $psiblast --cdhit $cdhit --blastdbcmd $blastdbcmd --num_threads 40 > ${i%.proveaninput.fa*}.proveanout"
  echo ${cmd}
  eval ${cmd}
  done
  cd
  
3. Sequence comparison between S.chilense and S.lycopersicoides
# The GO terms fasta file of S.chilense is the query file for the analysis
# Blast similarity search of GO terms in S.chilense against S.lycopersicoides
  $blast/makeblastdb -in $lycopersicoides -dbtype 'nucl'
  $blast/blastn -db $lycopersicoides -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cly/chilense.lycopersicoides.blastn

# Extract chilenseIDs and lycopersicoides IDs to a file
  awk '{print $1}' $cly/chilense.lycopersicoides.blastn > $cly/chilenseID_lycopersicoides
  awk '{print $2}' $cly/chilense.lycopersicoides.blastn > $cly/lycopersicoidesID

# Extract chilense and lycopersicoides IDS in each line to a separate file
  cd $cly
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.lycopersicoides.blastn

# Extract the aminoacid sequences
  for i in *.txt; 
  do 
  seqkit grep -f ${i} $wd/chilense_lycopersicoides_genome_aa.fa > ${i%.txt*}.aminoacid.fa; 
  done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
  for i in *.aminoacid.fa ; 
  do 
  awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; 
  done

# Store the S.chilense gene sequence from the above file to a separate file. (This is used for variant effect prediction analysis)  
  for i in *.new.fa; 
  do 
  awk '/^>/{if(N)exit;++N;} {print;}' ${i} > ${i%.new.fa*}.proveaninput.fa; 
  done
  
# Now perform MAFFT MSA
  for i in *.new.fa; 
  do 
  mafft ${i} > ${i%.new.fa*}.mafft_out.fa; 
  done
  
# Extract the positions of variants, insertions, and deletions from the MSA file. Remove thw whitespaces from file so that output is compatible for running PROVEAN. 
  for i in *.mafft_out.fa; 
  do 
  echo ${i} processing....
  python $wd/scripts/extract_variants_indels_from_MSA_v2.py ${i} > ${i%.mafft_out.*}.var
  done

# Run PROVEAN analysis
  for i in *.proveaninput.fa;
  do 
  cmd="$provean -q ${i} -d $nrdb -v ${i%.proveaninput.fa*}.var --psiblast $psiblast --cdhit $cdhit --blastdbcmd $blastdbcmd --num_threads 40 > ${i%.proveaninput.fa*}.proveanout"
  echo ${cmd}
  eval ${cmd}
  done
  cd 

4. Sequence comparison between S.chilense and S.pimpinelifolium (LA1589)
# The GO terms fasta file of S.chilense is the query file for the analysis
# Blast similarity search of GO terms in S.chilense against S.pimpinellifolium
  $blast/makeblastdb -in $pimpinellifolium -dbtype 'nucl'
  $blast/blastn -db $pimpinellifolium -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cpi/chilense.pimpinellifolium.blastn

# Extract chilenseIDs and pimpinellifolium IDs to a file
  awk '{print $1}' $cpi/chilense.pimpinellifolium.blastn > $cpi/chilenseID_pimpinellifolium
  awk '{print $2}' $cpi/chilense.pimpinellifolium.blastn > $cpi/pimpinellifoliumID

# Extract chilense and pimpinellifolium IDS in each line to a separate file
  cd $cpi
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.pimpinellifolium.blastn

# Extract the aminoacid sequences
  for i in *.txt; 
  do 
  seqkit grep -f ${i} $wd/chilense_pimpinellifolium_genome_aa.fa > ${i%.txt*}.aminoacid.fa; 
  done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
  for i in *.aminoacid.fa; 
  do 
  awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; 
  done

# Store the S.chilense gene sequence from the above file to a separate file. (This is used for variant effect prediction analysis)  
  for i in *.new.fa; 
  do 
  awk '/^>/{if(N)exit;++N;} {print;}' ${i} > ${i%.new.fa*}.proveaninput.fa; 
  done
  
# Now perform MAFFT MSA
  for i in *.new.fa; 
  do 
  mafft ${i} > ${i%.new.fa*}.mafft_out.fa; 
  done

# Extract the positions of variants, insertions, and deletions from the MSA file. Remove thw whitespaces from file so that output is compatible for running PROVEAN. 
  for i in *.mafft_out.fa; 
  do 
  echo ${i} processing....
  python $wd/scripts/extract_variants_indels_from_MSA_v2.py ${i} > ${i%.mafft_out.*}.var
  done

# Run PROVEAN analysis
  for i in *.proveaninput.fa;
  do 
  cmd="$provean -q ${i} -d $nrdb -v ${i%.proveaninput.fa*}.var --psiblast $psiblast --cdhit $cdhit --blastdbcmd $blastdbcmd --num_threads 40 > ${i%.proveaninput.fa*}.proveanout"
  echo ${cmd}
  eval ${cmd}
  done
  cd
  
5. Sequence comparison between S.chilense and S.sitiens
# The GO terms fasta file of S.chilense is the query file for the analysis
# Blast similarity search of GO terms in S.chilense against S.sitiens
  $blast/makeblastdb -in $sitiens -dbtype 'nucl'
  $blast/blastn -db $sitiens -query Results/chilense.GO.terms.salt.drought.IDs.fasta -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident evalue qcovs qcovhsp qlen slen qseq sseq stitle" > $cs/chilense.sitiens.blastn

# Extract chilenseIDs and sitiens IDs to a file
  awk '{print $1}' $cs/chilense.sitiens.blastn > $cs/chilenseID_sitiens
  awk '{print $2}' $cs/chilense.sitiens.blastn > $cs/sitiensID

# Extract chilense and sitiens IDS in each line to a separate file
  cd $cs
  awk '{print $1"\n"$2 > $1"_"$2 ".txt"}' chilense.sitiens.blastn

# Extract the aminoacid sequences
  for i in *.txt; 
  do 
  seqkit grep -f ${i} $wd/chilense_sitiens_genome_aa.fa > ${i%.txt*}.aminoacid.fa; 
  done

# The above command outputs multiple sequences in a single file. This is because, blastn outputs multiple hits for a single query sequence search. The multiple hits identified were duplicates. hence it needs to be removed, so that, each fasta file has no duplicate headers and no duplicate fasta sequences.
  for i in *.aminoacid.fa; 
  do 
  awk '/^>/{f=!d[$1];d[$1]=1}f' ${i} > ${i%.aminoacid.fa*}.new.fa; 
  done

# Store the S.chilense gene sequence from the above file to a separate file. (This is used for variant effect prediction analysis)  
  for i in *.new.fa; 
  do 
  awk '/^>/{if(N)exit;++N;} {print;}' ${i} > ${i%.new.fa*}.proveaninput.fa; 
  done
  
# Now perform MAFFT MSA
  for i in *.new.fa; 
  do 
  mafft ${i} > ${i%.new.fa*}.mafft_out.fa; 
  done

# Extract the positions of variants, insertions, and deletions from the MSA file. Remove thw whitespaces from file so that output is compatible for running PROVEAN. 
  for i in *.mafft_out.fa; 
  do 
  echo ${i} processing....
  python $wd/scripts/extract_variants_indels_from_MSA_v2.py ${i} > ${i%.mafft_out.*}.var
  done

# Run PROVEAN analysis
  for i in *.proveaninput.fa;
  do 
  cmd="$provean -q ${i} -d $nrdb -v ${i%.proveaninput.fa*}.var --psiblast $psiblast --cdhit $cdhit --blastdbcmd $blastdbcmd --num_threads 40 > ${i%.proveaninput.fa*}.proveanout"
  echo ${cmd}
  eval ${cmd}
  done
  cd
  
# Generate an Upset plot to get the gene IDs intersections between different species in the study. 
# A data table needs to be created with all the S.chilense IDs matched to the different species in the study. Once the data table is created, the R script to generate upset plot is saved as a separate script.
  paste $wd/Results/chilense_lycopersicum/chilenseID_lycopersicum $wd/Results/chilense_lycopersicoides/chilenseID_lycopersicoides $wd/Results/chilense_pennellii/chilenseID_pennellii $wd/Results/chilense_pimpinellifolium/chilenseID_pimpinellifolium $wd/Results/chilense_sitiens/chilenseID_sitiens > $wd/Results/GO.IDs.all.species

# Add column names
  sed  -i '1i chilense_lycopersicum\tchilense_lycopersicoides\tchilense_pennellii\tchilense_pimpinellifolium\tchilense_sitiens' $wd/Results/GO.IDs.all.species

# Run the Rscript to generate the upset plot
  Rscript $wd/scripts/create_upsetplot.R

# Generate a phylogenetic tree of all the Blast hits in lycopersicum, lycopersicoides, pennellii, pimpinellifolium, and sitiens against chilense GO terms
# Extract the fasta sequences of the IDs in the Blast hits
  seqkit grep -f $cl/lycopersicumID $lycopersicum | awk 'BEGIN {print ">S.lycopersicum"};!/^>/{print}' > $cl/lycopersicumID.fasta
  seqkit grep -f $cly/lycopersicoidesID $lycopersicoides | awk 'BEGIN {print ">S.lycopersicoides"};!/^>/{print}' > $cly/lycopersicoidesID.fasta
  seqkit grep -f $cp/pennelliiID $pennellii | awk 'BEGIN {print ">S.pennelli"};!/^>/{print}' > $cp/pennelliiID.fasta
  seqkit grep -f $cpi/pimpinellifoliumID $pimpinellifolium | awk 'BEGIN {print ">S.pimpinellifolium"};!/^>/{print}' > $cpi/pimpinellifoliumID.fasta
  seqkit grep -f $cs/sitiensID $sitiens | awk 'BEGIN {print ">S.sitiens"};!/^>/{print}' > $cs/sitiensID.fasta
  
 # Concatenate the above fasta sequences into a single multifasta file.
 # This is for running the Multiple Sequence Alignment
   cat $cl/lycopersicumID.fasta $cly/lycopersicoidesID.fasta $cp/pennelliiID.fasta $cpi/pimpinellifoliumID.fasta $cs/sitiensID.fasta > $wd/Results/LLPPS.fasta
 # count the number of sequences in the Multi fasta merged file. It should be same as the number of genomes in the study
   grep ">" $wd/Results/LLPPS.fasta
   
 # Run the Multiple Sequence Alignment program
   mafft $wd/Results/LLPPS.fasta > $wd/Results/LLPPS.mafft.out
   
 # Generate the Phylogenetic tree
   raxmlHPC -d -p 12345 -m GTRGAMMAI -s $wd/Results/LLPPS.mafft.out -n $wd/Results/LLPPS
   
 # Visualise the tree in any of the tree viewing softwares like Figtree or iTOL.
