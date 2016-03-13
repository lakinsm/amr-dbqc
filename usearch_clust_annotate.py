#!/usr/bin/env python3

## A quick script to annotate usearch -clust_fast -uc output with the database annotations.

#############
## Imports ##
#############
import argparse
import sys


##########
## Vars ##
##########
db_annot = {}


#############
## Methods ##
#############
def load_annots(annot_file):
    """
    Parses the production database annotation file and stores the
    gene headers as keys and the annotations as values in a list.
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
parser.add_argument('infile', type=str, help='File path to usearch cluster -uc output file')
parser.add_argument('annot_file', type=str, help='File path to the database annotation file')


##########
## Main ##
##########
args = parser.parse_args()
load_annots(args.annot_file)

with open(args.infile, 'r') as clstr:
    for line in clstr.read().strip().split('\n'):
        temp = line.split('\t')
        outline = '\t'.join([temp[0], temp[1], temp[2], temp[3], temp[8]] + db_annot[temp[8]])
        sys.stdout.write(outline + '\n')






