#!/usr/bin/env python3

#############
## Imports ##
#############
import argparse
import sys


##########
## Vars ##
##########
db_annot = {}
singletons = set()
pslx_hash = {}
db_hash = {}
count = 1
header_list = set()


#############
## Methods ##
#############
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


def pslx_parse(infile):
    """
    Parses a pslx file line by line.
    :param infile: path to the input pslx file
    :return: generator of relevant information
    """
    with open(infile, 'r') as pslx:
        while True:
            line = pslx.readline()
            if line.startswith('-'):
                break
        while True:
            line = pslx.readline()
            if not line:
                return  # Stop Iteration
            if line.startswith("#"):
                continue
            line = line.strip().split()
            data = (line[1],  # 0 mismatch
                    line[4],  # 1 query gap
                    line[6],  # 2 target gap
                    line[9],  # 3 query name
                    line[10],  # 4 query size
                    line[11],  # 5 query start
                    line[12],  # 6 query stop
                    line[13],  # 7 target name
                    line[14],  # 8 target size
                    line[15],  # 9 target start
                    line[16])  # 10 target stop
            yield data
        assert False, "Should not reach this line"


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
            if temp[0]:
                db_annot.setdefault(temp[0], temp[1:])


def load_singletons(sing_file):
    with open(sing_file, 'r') as sing:
        for line in sing:
            singletons.add(line.strip())


def determine_set(header1, header2):
    if header1 in singletons:
        return 1, header1
    else:
        return 2, header2


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('blat_truncate_singletons.py')
parser.add_argument('pslx_file', type=str, help='File path to pslx BLAT results')
parser.add_argument('database_annotations', type=str, help='File path to the database annotation csv')
parser.add_argument('database_seqs', type=str, help='File path to the database sequence fasta file')
parser.add_argument('singleton_headers', type=str, help='File path to singleton headers (no >)')
parser.add_argument('outannot', type=str, help='File path to output annotaiton file')
parser.add_argument('outseq', type=str, help='File path to output sequence file')


##########
## Main ##
##########
args = parser.parse_args()
load_annots(args.database_annotations)

for header, seq in fasta_parse(args.database_seqs):
    db_hash.setdefault(header, seq)

for entry in pslx_parse(args.pslx_file):
    qcov = float(int(entry[6]) - int(entry[5])) / int(entry[4])
    tcov = float(int(entry[10]) - int(entry[9])) / int(entry[8])
    if qcov >= 0.5 or tcov >= 0.5:
        truncate = determine_set(entry[3], entry[7])
        header, seq = (truncate[1], db_hash[truncate[1]])
        annot = db_annot[header]
        if truncate[0] == 1:
            seq = seq[int(entry[5]):int(entry[6])]
        else:
            seq = seq[int(entry[9]):int(entry[10])]
        header_list.add(header)
        newheader = header + '|singleton_truncated|' + str(count)
        if newheader in db_hash:
            newheader = newheader + '|' + str(count)
        db_annot[newheader] = annot
        db_hash[newheader] = seq
        count += 1

for name in header_list:
    del db_hash[name]
    del db_annot[name]

with open(args.outannot, 'w') as outa:
    for k, v in db_annot.items():
        outa.write(k + ',' + ','.join(v) + '\n')

with open(args.outseq, 'w') as outs:
    for k, v in db_hash.items():
        outs.write('>' + k + '\n' + v + '\n')


