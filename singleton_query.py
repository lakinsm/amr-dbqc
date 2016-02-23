#!/usr/bin/env python3


from Bio import SeqIO
from Bio import Entrez
import os

infile='biopython_search.txt'
Entrez.email = "Steven.Lakin@colostate.edu"     # Always tell NCBI who you are
counter = 1

with open(infile, 'r') as f, open('failed_biopython_queries.txt', 'w') as errfile,\
        open('entrez_genome_headers.txt', 'w') as genfile, open('entrez_parsed_seqs.fasta', 'w') as truseq,\
        open('entrez_parsed_metadata.txt', 'w') as truseq_meta:
    truseq_meta.write('number,query,accession,name,description,len,keywords')
    terms = f.read().strip().split()
    for i in terms:
        print(i)
        filename = i+'.gbk'
        if not os.path.isfile('entrez/'+filename):
            # Downloading...
            try:
                net_handle = Entrez.efetch(db="nucleotide",id=i,rettype="gb", retmode="text")
                out_handle = open('entrez/'+filename, "w")
                out_handle.write(net_handle.read())
                out_handle.close()
                net_handle.close()
                print("Saved")
            except Exception as e:
                print(e)
                errfile.write(i+'\n')
        else:
            print("Parsing...")
            record = SeqIO.read('entrez/'+filename, "genbank")
            if len(record.seq) > 20000:
                print('{:s} is a putative genome, skipping...\n'.format(i))
                genfile.write(i+'\n')
            else:
                for x in range(len(record.features)):
                    if record.features[x].type == 'CDS':
                        subset = record.features[x].extract(record)
                        truseq.write('>'+str(counter)+'|'+subset.id.replace('<unknown id>', 'unknown_id')+'|'+subset.name.replace('<unknown name>', 'unknown_name')+'\n'+str(subset.seq)+'\n')
                        truseq_meta.write(','.join([str(counter), i, subset.id.replace('<unknown id>', 'unknown_id'),
                                                    subset.name.replace('<unknown name>', 'unknown_name'),
                                                    subset.description.replace('<unknown description>', 'unknown_description'),
                                                    str(len(subset.seq)),
                                                    '|'.join(record.annotations['keywords'])])+'\n')
                        counter += 1
