#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Combines SL and Poly(A) site data across multiple stages to generate combined
files indicating the cumulative scores for each site.

Keith Hughitt (khughitt@umd.edu)
2016/04/11

Example Usage:
--------------

./combine_site_data.py tcruzi_infecting_hsapiens_amastigote_sl_sorted.gff.gz \
                       tcruzi_infecting_hsapiens_epimastigote_sl_sorted.gff.gz \
                       tcruzi_infecting_hsapiens_trypomastigote_sl_sorted.gff.gz

Note that the script currently assumes each of the input GFFs shares a common
filename prefix which is adapted for the output filename determination.

Source: https://github.com/elsayed-lab/utr_length_composition_and_site_usage_analysis
"""
import os
import sys
import gzip
from BCBio import GFF

def main():
    """Main script body"""
    # Do a quick check of inputs
    gffs = sys.argv[1:]

    for gff in gffs:
        if not os.path.isfile(gff):
            sys.exit("Invalid input file specified: %s" % gff)

    # Open first input GFF
    if gffs[0].endswith('.gz'):
        fp = gzip.open(gffs[0])
    else:
        fp = open(gffs[0])

    # Create a list to store output entries
    combined = []

    # Get GFF header and chromosome entries from first input file (it's the
    # same for all input files)
    for line in fp:
        if line.startswith("#") or "\tchromosome\t" in line:
            combined.append(line)
    
    # Reset read counter
    fp.seek(0)

    # Get entries for the first input file
    chromosomes = {}

    for entry in GFF.parse(fp):
        if len(entry.features) > 0 and entry.features[0].type in ['chromosome', 'contig']:
            chromosomes[entry.id] = entry

    fp.close()

    # Add sites from other input GFFs
    for gff in gffs[1:]:
        # Open first input GFF
        if gff.endswith('.gz'):
            fp = gzip.open(gff)
        else:
            fp = open(gff)

        for entry in GFF.parse(fp):
            for feature in entry.features:
                chromosomes[entry.id].features.append(feature)

    # combined sites
    sites = {}

    # Sort and combine sites and add to results    
    for ch_id in chromosomes:
        print("Parsing sites for %s" % ch_id)
        chromosomes[ch_id].features.sort()

        sites[ch_id] = {}

        for site in chromosomes[ch_id].features:
            # skip chromosomes
            if site.type in ['chromosome', 'contig']:
                continue

            # site info
            gene_id = site.qualifiers['Name'].pop()
            desc = ",".join(site.qualifiers['description'])
            score = int(site.qualifiers['score'].pop())
            source = site.qualifiers['source'].pop()
            feature_type = site.type

            if gene_id not in sites[ch_id]:
                sites[ch_id][gene_id] = {}
    
            # new entry            
            if site.location.end not in sites[ch_id][gene_id]:
                sites[ch_id][gene_id][site.location.end] = {
                    'Name': gene_id,
                    'source': source,
                    'type': feature_type,
                    'score': score,
                    'strand': site.strand,
                    'description': desc
                }
            else:
                # updating score
                sites[ch_id][gene_id][site.location.end]['score'] += score

    # Add combined rows
    for ch_id in sites:
        print("Combining sites for %s" % ch_id)
        # iterate over genes on each chromosome
        for gene_id in sites[ch_id]:
            site_counter = 1

            # iterate over sites for gene
            for loc in sorted(sites[ch_id][gene_id].keys()):
                site = sites[ch_id][gene_id][loc]

                # feature type abbreviation
                if site['type'] == 'trans_splice_site':
                    feature_type = 'sl'
                else:
                    feature_type = 'polya'
      
                # assign a new site id
                site_id = "%s.%s.%d" % (site['Name'], feature_type, site_counter)

                # description
                desc = "ID=%s;Name=%s;description=%s" % (site_id,
                                                         site['Name'],
                                                         site['description'])

                # create output row
                strand = '+' if site['strand'] == 1 else '-'

                row = "\t".join([ch_id, site['source'],
                                 site['type'], str(loc), str(loc),
                                 str(site['score']), strand, '.', desc])
                combined.append(row + "\n") 
                site_counter += 1

    # write combined gff
    #gff_suffix = commonsuffix(gffs).replace('.gz', '')
    gff_suffix = ".gff"

    outfile = "".join([os.path.commonprefix(gffs), 'combined', gff_suffix])

    with open(outfile, 'w') as output:
        output.writelines(combined)

    print("Done!")

def commonsuffix(strs):
    """Takes a list of strings and returns their common suffix"""
    reverse_strs = [x[::-1] for x in strs]
    return os.path.commonprefix(reverse_strs)[::-1]

if __name__ == '__main__':
    main()
