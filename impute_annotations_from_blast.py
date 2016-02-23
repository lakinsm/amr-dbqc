#!/usr/bin/env python3

#############
## Imports ##
#############
import argparse
import sys
from collections import Counter


##########
## Vars ##
##########
db_annot = {}
blast_hash = {}


#############
## Methods ##
#############
def blast_parse(infile):
    """
    Parses a pslx file line by line.
    :param infile: path to the input pslx file
    :return: generator of relevant information
    """
    with open(infile, 'r') as blast:
        while True:
            line = blast.readline()
            if not line:
                return  # Stop Iteration
            if line.startswith("#"):
                continue
            line = line.strip().split()
            yield line
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
            db_annot.setdefault(temp[0], temp[1:4])

##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('impute_annotations_from_blast.py')
parser.add_argument('blast_file', type=str, help='File path to format 6 (custom) BLAST results')
parser.add_argument('database_annotations', type=str, help='File path to the database annotation csv')


##########
## Main ##
##########
args = parser.parse_args()
load_annots(args.database_annotations)

# 0 query seq id
# 1 target accession
# 2 target seq id
# 3 query length
# 4 query start
# 5 query stop
# 6 target length
# 7 target start
# 8 target stop
# 9 evalue
# 10 percent identity
# 11 gaps
# 12 query coverage?
# 13 target all titles <> delim
# 14 target taxonomy ids
# 15 target common names
# 16 target blast names
# 17 target sequence

for entry in blast_parse(args.pslx_file):
    cov = float(int(entry[5]) - int(entry[4])) / int(entry[3])
    gaps = int(entry[11])
    ttitles = entry[13].replace('>', '|').replace('<', '').rstrip('|').split('|')
    ttax = entry[14].split(';')
    tcom = entry[15].split(';')
    tblast = entry[16].split(';')












