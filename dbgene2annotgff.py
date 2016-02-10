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
import sys


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
            if line.startswith('-'):
                break
        while True:
            line = pslx.readline()
            if not line:
                return  # Stop Iteration
            if line.startswith("#"):
                continue
            line = line.strip().split()
            # name, mismatch, qgap, tgap, strand, qname, qsize, qstart, qstop, tname, tstart, tstop, bsize
            data = (line[0],  # 0 name
                    line[1],  # 1 mismatch
                    line[4],  # 2 query gap
                    line[6],  # 3 target gap
                    line[8],  # 4 strand
                    line[9],  # 5 query name
                    line[10],  # 6 query size
                    line[11],  # 7 query start
                    line[12],  # 8 query stop
                    line[13],  # 9 target name
                    line[15],  # 10 target start
                    line[16],  # 11 target stop
                    line[18],  # 12 block size
                    line[19],  # 13 multi query starts
                    line[20],  # 14 multi target starts
                    line[21])  # 15 seqs
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
evals = [x for x in zip(*refseq_annot.values())][1]

sys.stdout.write('db_gene\tmismatch\tquery_gap\ttarget_gap\tstrand\tclass\tmechanism\tgroup\tfeature_start\tfeature_stop\tfeature_size\tfeature_string\tseqs\n')
for entry in pslx_parse(args.pslx_file):
    ## Start value
    annots = db_annot[entry[5]]
    hits = [x for x in zip(*[y.split(',') for y in entry[12:16]])]  # tuples of (size, qstart, tstart, seq) + empty
    for tup in hits:
        if all(x for x in tup):
            smatch = binary_search(svals, int(tup[2]))
            if smatch > int(tup[2]):
                sindex = svals.index(str(smatch)) - 1
            else:
                sindex = svals.index(str(smatch))
            ## End value
            tup_eval = int(tup[2]) + int(tup[0])
            ematch = binary_search(evals, tup_eval)
            if ematch < tup_eval:
                eindex = evals.index(str(ematch)) + 1
            else:
                eindex = evals.index(str(ematch))
            keyvalues = []
            for j in range(sindex, eindex + 1):
                keyvalues.append(next(v for i, v in enumerate(refseq_annot.items()) if i == j))
            for k, v in keyvalues:
                if int(v[0]) - int(tup[2]) <= 0:
                    covstart = int(tup[1])
                else:
                    covstart = int(v[0]) - int(tup[2])
                if int(v[1]) - tup_eval > 0:
                    covend = int(tup[1]) + int(tup[0])
                else:
                    covend = (int(tup[1]) + int(tup[0])) - (int(v[1]) - tup_eval)
                sys.stdout.write(str(entry[5])+'\t'+'\t'.join([str(x) for x in entry[2:5]])+'\t'+'\t'.join(annots)+'\t'
                                 +str(covstart)+'\t'+str(covend)
                                 +'\t'+str(tup[0])+'\t'+str(k)+'\t'+str(tup[3])+'\n')






