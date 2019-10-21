#parse-genbank.py modified from https://github.com/adina/scripts-for-ngs/blob/master/parse-genbank.py

import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


genome=SeqIO.read(sys.argv[1], 'genbank')


l = []
n = 0
for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
    name = record.name
    for feat in genome.features:
        if feat.type == "rRNA":
            if 'small subunit ribosomal RNA' in feat.qualifiers['product'][0]:#or '16S ribosomal' for strict match?
                start = feat.location.start.position
                end = feat.location.end.position
                pos = [start, end]
                l.append(pos)
                print '>' + name               

                print feat.extract(genome.seq)
                n = n + 1
