#!/usr/bin/env python3

## This script takes blast custom output and imputes a basic (high level) annotation based on the highest frequency
## of exact match in blast results.  This can optionally be paired with diamond blastx and blastdbcmd parse data
## for my custom annotation pipeline.  This script is highly specific for my use case.

#############
## Imports ##
#############
import argparse
import sys
import re
from collections import Counter


##########
## Vars ##
##########
db_annot = {}
diamond_hash = {}
db_hash = {}
diamond_reject = set()
annotated_seqs = set()
unannotated_seqs = set()


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
        progress = 0
        while True:
            progress += 1
            if progress % 200000 == 0:
                sys.stderr.write('\r{:d} hits processed'.format(progress))
            line = blast.readline()
            if not line:
                sys.stderr.write('\r{:d} hits processed\n'.format(progress))
                return  # Stop Iteration
            if line.startswith("#"):
                continue
            line = line.rstrip().split('\\t')
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


def load_diamond(diamond_file):
    """
    Parses diamond hit file (tabular blast output format) and stores the query
    :param diamond_file:
    :return:
    """
    with open(diamond_file, 'r') as dmnd:
        data = dmnd.read().split('\n')
        for line in data:
            temp = line.split('\t')
            if line:
                query = temp[1].split('|')
                if float(temp[2]) >= 70:
                    diamond_hash.setdefault((query[1], query[3]), [temp[0], Counter()])
                else:
                    diamond_reject.add((query[1], query[3]))


def pop_unknown(counter_obj):
    unknown = re.compile(r'unknown')
    org = re.compile(r'[\[](.*)[\]]$')
    try:
        annot = max(counter_obj)
        if unknown.search(max(counter_obj)):
            counter_obj.pop(annot)
            pop_unknown(counter_obj)
        else:
            try:
                origin = org.search(annot).group(1)
                annot = annot.replace('['+origin+']', '').rstrip()
            except AttributeError:
                origin = 'NA'
                annot = max(counter_obj)
            return annot, origin
    except ValueError:
        return None, None
    return None, None



##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('impute_annotations_from_blast.py')
parser.add_argument('diamond_file', type=str, help='File path to tab delimited DIAMOND blastx results')
parser.add_argument('blastdbcmd_file', type=str, help='File path to blastdbcmd file from diamond results')
parser.add_argument('database_annotations', type=str, help='File path to the database annotation csv')
parser.add_argument('unannotated_file', type=str, help='Output file to write remaining unannotated queries')


##########
## Main ##
##########
args = parser.parse_args()
load_annots(args.database_annotations)
load_diamond(args.diamond_file)
sys.stderr.write('{:d} accepted hits from diamond blastx\n'.format(len(diamond_hash)))
sys.stderr.write('{:d} rejected hits from diamond blastx\n'.format(len(diamond_reject)))

## Diamond blastx tabular output format
# 0 query ID
# 1 subject ID
# 2 percent identity
# 3 alignment length
# 4 mismatch count
# 5 gap open count
# 6 query start
# 7 query end
# 8 target start
# 9 target end
# 10 evalue
# 11 bit score

## Blastdbcmd file format (needs to be split on \\t)
# 0 sequence
# 1 protein identifier (e.g. WP_#########.#) accession
# 2 gene identifier (e.g. #########)
# 3 ordinal id (not really used)
# 4 full id
# 5 sequence title <--- important part for annotation of proteins using blastx
# 6 sequence length
# 7 taxonomy ID
# 8 node id
# 9 node id
# 10 parent node id
# 11 tax name
# 12 scientific tax name
# 13 blast tax name
# 14 parent tax name
# 15 PIG

## Annotation file format (comma-separated)
# 0 fasta header name
# 1 class
# 2 mechanism
# 3 group
# 4 PHScreen_Category
# 5 PHScreen_Type
# 6 PHScreen_Resistance

for entry in blast_parse(args.blastdbcmd_file):
    id = (entry[2], entry[1])
    if id in diamond_reject:
        continue
    elif id in diamond_hash:
        diamond_hash[id][1] += Counter((entry[5],))


for key, values in diamond_hash.items():
    if values[1]:
        annotated_seqs.add(values[0])

with open(args.unannotated_file, 'w') as un:
    for key, values in diamond_hash.items():
        annot, origin = pop_unknown(values[1])
        if annot:
            db_hash.setdefault(values[0], [Counter(), Counter(), Counter()])
            db_hash[values[0]][0] += Counter(('|'.join(key),))
            db_hash[values[0]][1] += Counter((annot,))
            db_hash[values[0]][2] += Counter((origin,))
        else:
            if values[0] in annotated_seqs:
                continue
            else:
                unannotated_seqs.add(values[0])
    for entry in unannotated_seqs:
        un.write(entry+'\n')
    for key, values in db_hash.items():
        try:
            sys.stdout.write(key+'\t'+'\t'.join([max(x) for x in values[1:]])+'\t'+'|'.join(values[1])+'\t'+';'.join(values[0])+'\n')
        except BrokenPipeError:
            break














