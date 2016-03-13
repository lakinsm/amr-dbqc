#!/usr/bin/env python3

## This script takes as input a database annotation file and removes duplicated header entries (these are
## accumulated over the course of annotations.  It will print to the screen if report is true and output to a new
## annotation file otherwise.


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
        data = annot.read().split('\n')[1:]
        for line in data:
            temp = line.split(',')
            if temp[0] and temp[0] != '':
                temp = ['' if x == 'NA' else x for x in temp]
                db_annot.setdefault(temp[0], []).append(temp[1:])


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('impute_annotations_from_blast.py')
parser.add_argument('infile', type=str, help='File path to dirty annotation file')
parser.add_argument('outfile', type=str, help='File path to output clean annotation file')
parser.add_argument('-r', '--report', action='store_true',
                    help='Output duplicated pairs to stdout, do not write to outfile')


##########
## Main ##
##########
args = parser.parse_args()
load_annots(args.infile)

with open(args.outfile, 'w') as out:
    for header, values in db_annot.items():
        if len(values) == 1:
            values[0][2] = values[0][2].upper()
            out.write(header+','+','.join(values[0])+'\n')
        else:
            if args.report:
                sys.stdout.write(header+'\t'+'\t'.join([','.join(x) for x in values])+'\n')
            else:
                temp = list(zip(*values))
                temp[2] = [x.upper() for x in temp[2]]
                values = list(zip(*temp))
                if len(set(values)) == 1:
                    out.write(header+','+','.join(values[0])+'\n')
                else:
                    l = min([x.count('') for x in values])
                    annot = [i for i, j in enumerate(values) if j.count('') == l]
                    if len(annot) > 1:
                        sys.stdout.write(header+'\t'+'\t'.join([','.join(values[x]) for x in annot])+'\n')
                        out.write(header+','+','.join(values[annot[1]])+'\n')  ## This line is specific for this run on March 5, 2016; delete if doing new cleaning
                    else:
                        out.write(header+','+','.join(values[annot[0]])+'\n')


