#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Performs a multiple sequence alignment on the predicted polypyrimidine tracts
for detected for each primary trans-splicing site.

Author: Keith Hughitt
Date: Aug 04, 2016

Usage: python polypyrimidine_msa.py path/to/polypyrimidine_tracts.csv
"""
import sys
import os
from csv import DictReader

def main():
    """Main application logic"""
    # check input filepath
    infile = sys.argv[1]
    if not os.path.exists(infile):
        print("Invalid filepath specified")
        sys.exit(1)

    # parse CSV
    ids = []
    seqs = []

    # generate multifasta from csv
    outfile = open(infile.replace('csv', 'fa'), 'wt')

    with open(infile) as fp:
        for entry in DictReader(fp):
            entry = '>{0}\n{1}\n'.format(entry['primary_sl'], entry['seq'])
            outfile.write(entry)

    # run MAFFT

if __name__ == '__main__':
    main()


