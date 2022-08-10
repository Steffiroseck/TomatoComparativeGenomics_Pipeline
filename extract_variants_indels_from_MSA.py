import sys
from Bio import SeqIO

f = sys.argv[1]

# reading the query sequence (lycopersicum protein sequence line by line)
seq_records = SeqIO.parse(f, 'fasta')
# reading the reference sequence (chilense sequence line by line)
refseq_record = next(seq_records)
seq_record = next(seq_records)
# parsing the lycopersicum sequence first and checking if a.a are substituted or inserted or deleted.

index_iter = iter(range(0,len(refseq_record)))
ref_real_pos = 0
query_real_pos =0
for i in index_iter:
    # each amino acids in the chilense are stored one per line in nt1 and lycopersicum are stored in nt2 variables.
    nt1 = refseq_record[i]
    nt2 = seq_record[i]
    if nt1 != "-" and nt2 != "-":
        ref_real_pos = ref_real_pos + 1
        query_real_pos = query_real_pos + 1
    if nt1 == "-":
        deletion = ""
        while refseq_record[i] == "-":
            next (index_iter)
            nt1 = refseq_record[i]
            nt2 = seq_record[i]
            deletion = deletion + nt2 + str(query_real_pos + 1) + "_" 
            i = i+1
            query_real_pos = query_real_pos + 1
        print(deletion[0:len(deletion)-1]+"del")
    if nt2 == "-":
        insertion = ""
        start= i
        ref_real_start = ref_real_pos
        while seq_record[i] == "-":
            next (index_iter)
            nt1 = refseq_record[i]
            #nt2 = seq_record[i]
            insertion = insertion + nt1 
            i = i+1
            ref_real_pos = ref_real_pos + 1
        end = i
        ref_real_end = ref_real_pos
        print(seq_record[start-1]+str(ref_real_start)+"_"+seq_record[end]+str(ref_real_end+1)+"ins"+insertion)
    if nt1 != nt2 and nt1 != "-" and nt2 != "-":
        print(nt1+ str(ref_real_pos)+ nt2)
