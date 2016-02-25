#!/usr/bin/env python3

## First pass at annotation for post_diamond majority vote annotation file
## see the final_qc/blast_parsing folder for details

mechfile = 'post_diamond_mechlist.txt'
seq_annotfile = 'post_diamond_annotated_seqs.tsv'
outfile = 'post_diamond_passes/ready_for_review_annotations_pass1.csv'
leftfile = 'post_diamond_passes/remaining_seqs_pass1.tsv'

mech_hash = {}

with open(mechfile, 'r') as mech, open(seq_annotfile, 'r') as seq, open(outfile, 'w') as out, open(leftfile, 'w') as left:
    for line in mech:
        if line:
            data = line.strip().split('\t')
            mech_hash.setdefault(data[1], data[0])
    for line in seq:
        if line:
            data = line.strip().split('\t')
            guess = mech_hash[data[1]]
            if guess == 'NaN' or guess == 'NRes':
                left.write('\t'.join(data)+'\n')
            else:
                out.write(data[0]+','+guess+',,,\n')



