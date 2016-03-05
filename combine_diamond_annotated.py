#!/usr/bin/env python3

## First pass at annotation for post_diamond majority vote annotation file
## see the final_qc/blast_parsing folder for details

infile = 'ready_for_review_annotations_pass1.csv'
inseq = '/home/lakinsm/Documents/phd/phdenv/final_model_qc/todo_files/blast_parsing/data/fixed_full_blast_query.fasta'
inannot = 'hmm_annotations3.csv'
outseq = 'temp_hmm_seqs4.fasta'
outannot = 'temp_hmm_annotations4.csv'

db_annot = {}
db_seqs = {}

def load_annots(annot_file):
    """
    Parses the production database annotation file and stores the
    gene headers as keys and the annotations as values.
    :param annot_file: production database annotation file
    :return: void
    """
    with open(annot_file, 'r') as annot:
        data = annot.read().split('\n')
        for line in data:
            temp = line.split(',')
            db_annot.setdefault(temp[0], temp[1:])


def fasta_parse(infile):
    """ Parses a fasta file in chunks of 2 lines.
    :param infile: path to the input fasta file
    :return: generator of (header, sequence) fasta tuples
    """
    with open(infile, 'r') as fasta_file:
        # Skip whitespace
        while True:
            line = fasta_file.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            all_lines = []
            line = fasta_file.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                all_lines.append(line.rstrip())
                line = fasta_file.readline()
            yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"


load_annots(inannot)
with open(infile, 'r') as inf, open(outseq, 'w') as outs, open(outannot, 'w') as outa:
    for line in fasta_parse(inseq):
        db_seqs.setdefault(line[0], line[1])
    for line in inf:
        line = line.strip().split(',')
        if line[0] in db_annot:
            if line[1:].count('') >= db_annot[line[0]].count(''):
                annot = db_annot[line[0]]
            else:
                annot = line[1:]
                annot += [''] * (6 - len(annot))
        else:
            annot = line[1:]
            annot += [''] * (6 - len(annot))
        outs.write('>'+line[0]+'\n'+db_seqs[line[0]]+'\n')
        outa.write(line[0]+','+','.join(annot)+'\n')
