#!/usr/bin/env python3

## This script attempts to impute annotations by majority vote; if there is a choice between two, it prints to
## screen and doesn't annotate that level


import re
import argparse


db_annot = {}


parser = argparse.ArgumentParser('cdhit_impute_annotations.py')
parser.add_argument('clstr_file', type=str, help='File path to output from cdhit')
parser.add_argument('annotation_file', type=str, help='File path to MEG AMR annotation file')
parser.add_argument('outfile', type=str, help='File path to output annotation file')
parser.add_argument('-v', '--verbose', action='store_true', help='Print newline separators between clusters in outfile')


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


def separate_clusters(infile, outfile):
    cluster_num = 0
    clstr = []
    annot = {}
    header_reg = re.compile(r'>(.+?)\.{3}')
    with open(infile, 'r') as f, open(outfile, 'w') as out:
        line = f.readline()
        while line:
            if line[0] is ">":
                cluster_num += 1
                if cluster_num > 1:
                    if len(clstr) is 1:
                        for header in clstr:
                            out.write(header + ',' + ','.join([str(x) for x in db_annot[header]]) + '\n')
                            db_annot.pop(header, None)
                            if args.verbose:
                                out.write('\n')
                        clstr = []
                    else:
                        for header in clstr:
                            annot.setdefault(header, db_annot[header])
                        for i, level in enumerate(zip(*annot.values())):
                            if i > 2:
                                continue
                            val = set([x for x in level if x])
                            if len(val) == 1:
                                for k, v in annot.items():
                                    if not annot[k][i]:
                                        annot[k][i] = list(val)[0]
                            elif args.verbose:
                                print(list(annot.keys()), val)
                        for header in clstr:
                            out.write(header + ',' + ','.join([str(x) for x in annot[header]]) + '\n')
                            db_annot.pop(header, None)
                        if args.verbose:
                            out.write('\n')
                        clstr = []
                        annot = {}
            else:
                clstr.append(header_reg.findall(line)[0])
            line = f.readline()
        for k, v in db_annot.items():
            out.write(k + ',' + ','.join(v) + '\n')


args = parser.parse_args()
load_annots(args.annotation_file)

separate_clusters(args.clstr_file, args.outfile)
