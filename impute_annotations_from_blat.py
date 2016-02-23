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
pslx_hash = {}


#############
## Methods ##
#############
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
            if int(line[5]) > 100 or int(line[7]) > 100:
                continue
            data = (line[1],  # 0 mismatch
                    line[4],  # 1 query gap
                    line[6],  # 2 target gap
                    line[9],  # 3 query name
                    line[10],  # 4 query size
                    line[11],  # 5 query start
                    line[12],  # 6 query stop
                    line[13])  # 7 target name
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
            db_annot.setdefault(temp[0], temp[1:4])

##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('impute_annotations_from_blat.py')
parser.add_argument('pslx_file', type=str, help='File path to pslx BLAT results')
parser.add_argument('database_annotations', type=str, help='File path to the database annotation csv')


##########
## Main ##
##########
args = parser.parse_args()
load_annots(args.database_annotations)

for entry in pslx_parse(args.pslx_file):
    coverage = float(int(entry[6]) - int(entry[5])) / int(entry[4])
    mis = float(entry[0]) / int(entry[4])
    gap = float(entry[1]) / int(entry[4])
    pslx_hash.setdefault(entry[3].split('|')[0], []).append((coverage, mis, gap, entry[7], entry[3]))

for key, value in pslx_hash.items():
    annot_list = []
    mval = max(tuple(zip(*value))[0])
    ind = [i for i, v in enumerate(tuple(zip(*value))[0]) if v == mval]
    for x in ind:
        annot_list.append(tuple(db_annot[value[x][3]]))
    cannot = Counter(annot_list)
    sys.stdout.write(value[0][4]+','+','.join([str(x) for x in max(cannot)])+'\n')








