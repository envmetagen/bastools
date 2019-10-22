#parse-genbank.py modified from https://github.com/adina/scripts-for-ngs/blob/master/parse-genbank.py

import sys 
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


genome=SeqIO.read(sys.argv[1], 'genbank')


l = []
n = 0
for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
    name = record.name
    for feature in record.features:
      if 'source' in feature.type:
        taxid=''.join(feature.qualifiers['db_xref'])
        taxid=re.sub(r'.*taxon:','',taxid)
        org=''.join(feature.qualifiers['organism'])
        org=re.sub(r'.*organism=','',org)
    for feat in genome.features:
        if feat.type == "CDS":
            if 'cytochrome oxidase subunit 1' in feat.qualifiers['product'][0]: 
            # POSSIBLE VALUES FOR if
              #cytochrome oxidase subunit 1
              #cytochrome c oxidase subunit 1
              #cytochrome c oxidase subunit I
              #cytochrome oxidase subunit I
            #POSSIBLE VALUES FOR feat.qualifiers
              #product
              #note
                start = feat.location.start.position
                end = feat.location.end.position
                pos = [start, end]
                l.append(pos)
                print '>' + name + ' organism=' + org + '; taxid=' + taxid + ";"                 

                print feat.extract(genome.seq)
                n = n + 1
