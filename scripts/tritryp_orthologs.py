#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Determines three way orthologs for a set of three organisms in TriTrypDB
Keith Hughitt
2016/12/15
"""
import os
import gzip

# TriTrypDB Gene txt file for one of the three organisms
txt = os.path.join(os.getenv('REF'),
        'tcruzi_clbrener_esmeraldo-like/annotation',
        'TriTrypDB-28_TcruziCLBrenerEsmeraldo-likeGene.txt.gz')

# Two other organisms to check for orthologs with
org1 = 'Leishmania major strain Friedlin'
org2 = 'Trypanosoma brucei TREU927'

# iterate over file and parse ortholog tables
gene_id = None
section = 'main'

# list to keep track of three-way orthologs
orthologs = []

# keep track of number of orthologs
org1_ortholog = False
org2_ortholog = False

fp = gzip.open(txt, 'rt')
for line in fp:
    if section == 'main':
        # parse gene ids
        if line.startswith('Gene ID'):
            gene_id = line.split(':').pop().strip()
            print("Parsing gene %s" % gene_id)

        # start parsing ortholog table
        elif line.startswith('TABLE: Orthologs'):
            section = 'orthologs'
    elif section == 'orthologs':
        # parse ortholog table
        # skip header column
        line = next(fp)

        # check for both other organism
        while line != '\n':
            organism = line.split('\t')[1] 

            if organism == org1:
                org1_ortholog = True
            elif organism == org2:
                org2_ortholog = True

            line = next(fp)

        # if both are found, store gene id
        if org1_ortholog and org2_ortholog:
            orthologs.append(gene_id)

        # done parsing table
        org1_ortholog = False
        org2_ortholog = False
        section = 'main'

fp.close()

# write results
with open('tritryp_orthologs.txt', 'w') as fp:
    fp.write('\n'.join(orthologs))


