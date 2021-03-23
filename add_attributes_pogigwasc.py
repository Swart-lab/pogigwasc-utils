#!/usr/bin/env python

import re
import json
import argparse

from collections import defaultdict
from pogigwasc_utils_shared import *

def create_parent_genes(features):
    # Create parent genes object and change features object in place
    geneid_counter = defaultdict(int)
    genes = defaultdict( #seqid
        lambda: defaultdict( # geneid
            dict
        ))
    for seqid in features:
        features_list = [(f, features[seqid][f]) for f in features[seqid] if
                features[seqid][f]['type'] in ['CDS','stop_codon','intron']]
        features_list = sorted(features_list, key=lambda f: f[1]['start'])
        feature_ids_sorted = [f[0] for f in features_list]
        last_cluster = []

        for i in range(1, len(feature_ids_sorted)):
            curr_feature = feature_ids_sorted[i]
            last_feature = feature_ids_sorted[i-1]
            if features[seqid][last_feature]['end'] == features[seqid][curr_feature]['start'] \
                and features[seqid][last_feature]['strand'] == features[seqid][curr_feature]['strand']:
                feature_strand = features[seqid][curr_feature]['strand']
                # Adjacent features belong to same gene
                # Check for CDS adjacent to stop codon, adjust CDS coords 
                # to encompass stop codon
                if features_list[i][1]['type'] == 'CDS' \
                    and features_list[i-1][1]['type'] == 'stop_codon' \
                    and feature_strand == '-' :
                    # stop codon followed by CDS
                    features_list[i][1]['start'] = features_list[i-1][1]['start']
                elif features_list[i][1]['type'] == 'stop_codon' \
                    and features_list[i-1][1]['type'] == 'CDS' \
                    and feature_strand == '+' :
                    # CDS followed by stop codon
                    features_list[i-1][1]['end'] = features_list[i][1]['end']
                last_cluster.append(last_feature)
            else:
                # close current cluster and start a new one
                last_cluster.append(last_feature)
                geneid_counter[seqid] += 1
                geneid = seqid + ".g" + str(geneid_counter[seqid])
                genes[seqid][geneid] = {f:1 for f in last_cluster}
                last_cluster = []
        # close the last cluster
        last_cluster.append(curr_feature)
        geneid_counter[seqid] += 1
        geneid = seqid + ".g" + str(geneid_counter[seqid])
        genes[seqid][geneid] = {f:1 for f in last_cluster}
    return(genes)


if __name__ == "__main__":
    # command line args
    parser = argparse.ArgumentParser(
        description="""
        Create parent gene objects and attributes for Pogigwasc GFF output.

        Add attributes (ID, Parent), creates gene objects for adjacent
        CDS/intron/stop_codon features, adjust CDS coordinates to include
        adjacent stop codons, and add phase field to GFF for CDS features.
        """)
    parser.add_argument("-i", "--input", action='store', help="""
        Input GFF file from Pogigwasc
        """)
    parser.add_argument("-o", "--output", action='store', help="""
        Output GFF file
        """)
    args = parser.parse_args()

    # main script
    features = read_pogigwasc_gff(args.input)
    genes = create_parent_genes(features)
    out = []
    for seqid in genes:
        for geneid in genes[seqid]:
            gene_start, gene_end, gene_strand = find_gene_limits(genes, seqid, geneid, features)
            gene_start, gene_end = py2gff_coords(gene_start, gene_end)
            out.append([seqid, 'predicted','gene',gene_start, gene_end, '.', gene_strand, '.', f'ID={geneid}'])
            for feature_id in genes[seqid][geneid]:
                new_feature_id = feature_id
                # Renumber CDSs, because multiple CDS segments belonging to the
                # same gene should all have the same CDS ID although we earlier
                # used distinct identifiers to have a unique handle
                if features[seqid][feature_id]['type'] == 'CDS':
                    new_feature_id = seqid + '.c' + geneid.split(".")[-1][1:]
                out.extend(dict2gff(features[seqid][feature_id], new_feature_id, geneid))
    # write new file
    with open(args.output, 'w') as fh:
        for line in out:
            fh.write("\t".join([str(i) for i in line]))
            fh.write("\n")
