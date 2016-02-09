#!/usr/bin/env python3

## This script takes as input the annotations for the MEG AMR database, a refseq annotation file in gff3 format,
## and the results from a blat against the corresponding refseq genome in pslx format.  It outputs analytic data
## in csv format including all relevant annotation information as well as % feature coverage relative to the gene
## of interest from the AMR database.  This script was used for database quality control in order to find
## AMR homologue sequences that are problematic due to point mutations.

#############
## Imports ##
#############
import argparse
import collections


##########
## Vars ##
##########
db_annot = {}
refseq_annot = collections.OrderedDict()
memory = 0
pval = None


#############
## Methods ##
#############
def pslx_parse(infile):
    with open(infile, 'r') as pslx:
        """
        Parses a pslx file line by line.
        :param infile: path to the input pslx file
        :return: generator of relevant information
        """
        while True:
            line = pslx.readline()
            if line.startswith("#"):
                continue
            line = line.strip().split()
            # name, mismatch, qgap, tgap, strand, qstart, qsize, qstop, tname, tstart, tstop, bsize
            data = (line[0], line[1], line[4], line[6], line[8], line[10], line[11], line[12],
                    line[13], line[14], line[15], line[16], line[18])
            yield data
            if not line:
                return  # Stop Iteration
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
            db_annot.setdefault(temp[0], temp[1:])


def load_gff(gff_file):
    """
    Parses a gff file and stores the necessary information for further
    annotation of the AMR db genes.
    :param gff_file: target gff3 file
    :return: void
    """
    with open(gff_file, 'r') as gff:
        data = gff.read().split('\n')
        for line in data:
            if not line or line[0] == '#':
                continue
            temp = line.split('\t')
            if temp[1] == "ena" and temp[2] == "gene":
                refseq_annot.setdefault(temp[8], (temp[3], temp[4]))


def binary_search(val_list, term):
    """
    Standard divide and conquer binary search algorithm.
    Returns the match or one of the two neighbors.
    :param val_list: list of numeric integers to search
    :param term: numeric query value (integer)
    :return: match or value of one of the two neighbors
    """
    if len(val_list) == 1:
        return int(val_list[0])
    else:
        midpoint = len(val_list) // 2
        if int(val_list[midpoint]) == term:
            return int(val_list[midpoint])
        else:
            if term < int(val_list[midpoint]):
                return binary_search(val_list[:midpoint], term)
            else:
                return binary_search(val_list[midpoint + 1:], term)



##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('dbgene2annotgff.py.py')
parser.add_argument('pslx_file', type=str, help='File path to pslx BLAT results')
parser.add_argument('database_annotations', type=str, help='File path to the database annotation csv')
parser.add_argument('gff_annotations', type=str, help='File path to the GFF annotation file')


##########
## Main ##
##########
args = parser.parse_args()
load_gff(args.gff_annotations)
load_annots(args.database_annotations)
svals = [x for x in zip(*refseq_annot.values())][0]

for entry in pslx_parse(args.pslx_file):
    match = binary_search(svals, entry[10])
    if match > int(entry[10]):
        index = svals.index(str(match)) - 1
    else:
        index = svals.index(str(match))
        key, value = next(v for i, v in enumerate(refseq_annot.items()) if i == index)
    ## For reference in entry
    # 0 name
    # 1 mismatch
    # 2 qgap
    # 3 tgap
    # 4 strand
    # 5 qstart
    # 6 qsize
    # 7 qstop
    # 8 tname
    # 9 tstart
    # 10 tstop
    # 11 bsize




