#!/usr/bin/env python3

import re
from collections import defaultdict
import argparse
import json

from pogigwasc_utils_shared import *

def update_feature_coords_realtrons(features, introns):
    # Adjust coordinates of features with introns added
    # Record features that contain introns
    # features object is changed in-place

    intron_parents = defaultdict( # seqid
        lambda: defaultdict( # feature id
            list))

    for seqid in introns:
        try:
            # sort features ascending by start coordinate
            features_list = [ (f, features[seqid][f]) for f in features[seqid] ]
            features_list = sorted(features_list, key=lambda f: f[1]['start'])
            feature_ids_sorted = [f[0] for f in features_list]

            # sort introns ascending by start coordinate
            introns_list = [ (i, introns[seqid][i]) for i in introns[seqid]]
            introns_list = sorted(introns_list, key=lambda i: i[1]['start'])
            intron_ids_sorted = [i[0] for i in introns_list]

            # Adjust coordinates
            for intron_id in intron_ids_sorted:
                intron_len = introns[seqid][intron_id]['end'] - introns[seqid][intron_id]['start']

                for feature_id in feature_ids_sorted:
                    # sanity check
                    orig_feature_length = features[seqid][feature_id]['end'] - features[seqid][feature_id]['start']
                    # sequentially iterate through features and update coordinates
                    # assume that no two genes share the same start coordinate
                    if features[seqid][feature_id]['start'] > introns[seqid][intron_id]['start']:
                        # feature to right of intron: push whole feature forward
                        features[seqid][feature_id]['start'] += intron_len
                        features[seqid][feature_id]['end'] += intron_len
                    # in edge case where feature start == intron start, the
                    # resulting feature will have correct length, but the 
                    # first feature segment will have zero length
                    # however we want to record the corresponding gene as parent
                    # to the intron
                    else:
                        if features[seqid][feature_id]['end'] >= introns[seqid][intron_id]['start']:
                            # feature start to left of intron, feature start to right of intron
                            # i.e. intron contained in feature, need to adjust:
                            # keep start coordinate, push end coordinate forward
                            features[seqid][feature_id]['end'] += intron_len
                            intron_parents[seqid][feature_id].append(intron_id)
                            # sanity check
                            if features[seqid][feature_id]['end'] - features[seqid][feature_id]['start'] != orig_feature_length + intron_len:
                                print("feature length does not match after adding introns")
                        else:
                            # entire feature to left of intron
                            pass

        except KeyError:
            print(f"seqid {seqid} in intron gff not found in feature gff")
            pass
    return(intron_parents)


def check_strand_conflict_orphan_introns(features, introns, intron_parents):
    # Report adjusted coordinates
    # TODO: what about features that encompass both conflicting and non-conflicting introns?

    orphan_introns = defaultdict(list)
    conflict_features = defaultdict( # seqid
        lambda: defaultdict( # feature id
            list))

    for seqid in features:
        # sort introns ascending by start coordinate
        introns_list = [ (i, introns[seqid][i]) for i in introns[seqid] ]
        introns_list = sorted(introns_list, key=lambda i: i[1]['start'])
        intron_ids_sorted = [i[0] for i in introns_list]

        # Check for mismatched strands in intron vs. CDS
        for feature_id in intron_parents[seqid]:
            parent_strand = features[seqid][feature_id]['strand']
            for intron_id in intron_parents[seqid][feature_id]:
                if introns[seqid][intron_id]['strand'] != parent_strand:
                    # feature is a 'conflicted' CDS
                    conflict_features[seqid][feature_id].append(intron_id)

        # Check for intron orphans
        contained_introns = []
        for i in intron_parents[seqid].values():
            contained_introns.extend(i)
        orphan_introns[seqid].extend([i for i in intron_ids_sorted if i not in set(contained_introns)])
    return(orphan_introns, conflict_features)


def create_parent_genes(features):
    # Associate stop_codons with adjacent attached CDSs in parent gene objects
    # add attributes field to the feature table
    # assumes that:
    #  - every CDS is adjacent to a stop codon,
    #  - strand of CDS and adjacent stop codon always the same
    #  - that a gene is composed of no more than one CDS+stop codon

    geneid_counter = defaultdict(int)
    genes = defaultdict( # seqid
        lambda: defaultdict( # geneid
            dict)) # dict of features keyed by ID

    for seqid in features:
        # sort features by start coord
        features_list = [ (f, features[seqid][f]) for f in features[seqid] ]
        features_list = sorted(features_list, key=lambda f: f[1]['start'])
        feature_ids_sorted = [f[0] for f in features_list]

        last_feature = None # record previous feature

        for feature_id in feature_ids_sorted:
            if last_feature:
                # adjacent features
                if int(features[seqid][last_feature]['end']) == int(features[seqid][feature_id]['start']):
                    # either: CDS followed by stop codon and strands both positive
                    # or: stop codon followed by CDS, and strands both negative
                    if (features[seqid][last_feature]['type'] == 'CDS' \
                            and features[seqid][feature_id]['type'] == 'stop_codon' \
                            and features[seqid][last_feature]['strand'] == '+' \
                            and features[seqid][feature_id]['strand'] == '+'):
                        # extend CDS to encompass stop codon limits
                        features[seqid][last_feature]['end'] = features[seqid][feature_id]['end']
                        # update gene counter
                        geneid_counter[seqid] += 1
                        geneid = "_".join([seqid, 'gene', str(geneid_counter[seqid])])
                        # add feature ID and Parent gene id to attributes field
                        genes[seqid][geneid] = { feature_id : 1,
                                                 last_feature: 1}
                    elif (features[seqid][last_feature]['type'] == 'stop_codon' \
                        and features[seqid][feature_id]['type'] == 'CDS' \
                        and features[seqid][last_feature]['strand'] == '-' \
                        and features[seqid][feature_id]['strand'] == '-'):
                        # extend CDS to encompass stop codon limits
                        features[seqid][feature_id]['start'] = features[seqid][last_feature]['start']
                        # update gene counter
                        geneid_counter[seqid] += 1
                        geneid = "_".join([seqid, 'gene', str(geneid_counter[seqid])])
                        # add feature ID and Parent gene id to attributes field
                        genes[seqid][geneid] = { feature_id : 1,
                                                 last_feature: 1}
                    else:
                        print("adjacent features but apparently unrelated")
                # Update previous feature
                last_feature = feature_id
            else:
                # first entry
                last_feature = feature_id
    return(genes)


def segment_cds_with_introns(genes, features, introns, intron_parents, orphan_introns, conflict_features):
    # Changes `features` and `introns` in place

    for seqid in genes:
        # Split CDS features into segments and correct phase
        #  pogigwasc in intronless mode produces CDS features with phase 0 (check!)
        #  so the first segment will always have phase 0
        # Report features in GFF format
        for geneid in genes[seqid]:
            for feature_id in genes[seqid][geneid]:
                # Features that contain introns
                if feature_id in intron_parents[seqid]:

                    # strandedness of CDS and intron features conflict
                    if feature_id in conflict_features[seqid]:
                        # Note that strand is conflicting
                        features[seqid][feature_id]['Note'] = 'strand_conflict'
                        # report overlapping introns as 'intron_orphan' features
                        for conflict_intron in conflict_features[seqid][feature_id]:
                            introns[seqid][conflict_intron]['Note'] = 'orphan'

                    else:
                        # check if CDS feature (or start_codon or stop_codon) and split if necessary
                        if features[seqid][feature_id]['type'] in ['CDS', 'stop_codon']:
                            # get coordinates of the split CDS feature segments,
                            # interrupted by introns
                            seg_coords = [features[seqid][feature_id]['start']]
                            for intron_id in intron_parents[seqid][feature_id]:
                                # the introns should be present in coordinate sort order
                                seg_coords.append(introns[seqid][intron_id]['start'])
                                seg_coords.append(introns[seqid][intron_id]['end'])
                                # add gene as parent to this intron
                                introns[seqid][intron_id]['Parent'] = geneid
                            seg_coords.append(features[seqid][feature_id]['end'])

                            # start and end coords as a list for each segment
                            starts = seg_coords[::2]
                            ends = seg_coords[1::2]
                            segs = []
                            # forward strand
                            if features[seqid][feature_id]['strand'] == '+':
                                segs = list(zip(starts,ends))
                            # reverse strand: work backwards
                            elif features[seqid][feature_id]['strand'] == '-':
                                segs = reversed(list(zip(starts,ends)))
                            # replace feature 'start', 'end', and 'phase' fields with lists
                            features[seqid][feature_id]['start'] = starts
                            features[seqid][feature_id]['end'] = ends
                            features[seqid][feature_id]['phase'] = [0] * len(starts)
                        else:
                            print(f"Non CDS feature {feature_id} overlaps with empirical intron")
                            # this should never happen because update_feature_coords_realtrons() 
                            # is called with pogigwasc feature table which has only CDS and stop_codons
                else:
                    # intronless feature, report as-is, change phase from '.' to 0 if CDS
                    if features[seqid][feature_id]['type'] == 'CDS':
                        features[seqid][feature_id]['phase'] = 0
        # add orphan introns
        for orphan_id in orphan_introns[seqid]:
            introns[seqid][orphan_id]['Note'] = 'orphan'


def dict2gff(feature: dict, feature_id, parent_id=None):
    """Convert dict of features to lists for GFF output

    If feature is multisegmented, i.e. fields 'start', 'end', and 'phase' are
    list objects, then report each segment as separate line. 

    Parameters
    ----------
    feature : dict
        Features keyed by GFF keys
    feature_id : str
        Feature ID to use in ID field of attributes
    parent_id : str
        Parent feature ID to use in Parent field of attributes

    Returns
    -------
    list
        list of lists, containing 9 fields of GFF line
    """
    # return list of lists of GFF fields
    out = []
    #attrs = {'ID': feature_id, 'Parent' : parent_id}
    attrs = defaultdict(str)
    if 'attributes' in feature:
        for i in feature['attributes'].rstrip(';').split(';'):
            attrs[i.split('=')[0]] = attrs[i.split('=')[1]]
    attrs['ID'] = feature_id
    if parent_id:
        attrs['Parent'] = parent_id
    elif 'Parent' in feature:
        attrs['Parent'] = feature['Parent']
    # add other attributes
    attrkeys = ['ID','Parent']
    for key in feature:
        if key not in GFF_KEYS and key not in ['ID','Parent']:
            attrs[key] = feature[key]
            attrkeys.append(key) # to ensure that ID and Parent come first
    attributes_string = ";".join([str(key)+"="+str(attrs[key]) for key in attrkeys if key in attrs])

    if type(feature['start']) is list:
        # multi-segment feature, all fields identical except start/end/phase coords
        for i in range(len(feature['start'])):
            # skip zero-length feature segments
            if feature['end'][i] > feature['start'][i]:
                newstart, newend = py2gff_coords(feature['start'][i], feature['end'][i])
                line = [feature['seqid'], feature['source'], feature['type'],
                        newstart, newend,
                        feature['score'], feature['strand'],
                        feature['phase'][i],
                        attributes_string]
                out.append(line)
            elif feature['start'][i] == feature['end'][i]:
                print(f"Skipping zero-length feature segment { str(feature['start'][i]) } for feature {feature_id}")
            else:
                print(f"Start before end for coords { str(feature['start'][i]) } { str(feature['end'][i]) } for feature {feature_id}")

    else:
        newstart, newend = py2gff_coords(feature['start'], feature['end'])
        line = [feature['seqid'], feature['source'], feature['type'],
                newstart, newend,
                feature['score'], feature['strand'],
                feature['phase'],
                attributes_string]
        out.append(line)

    return(out)


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--introns", action='store', type=str, 
        help="Realtrons GFF3 file")
    parser.add_argument("--features", action='store', type=str, 
        help="""
        Pogigwasc intronless-mode gene predictions, GFF3 file.
        Assumes no overlapping features, only `CDS` and `stop_codon` features,
        and stop codons not contained in adjacent CDS features.
        """)
    parser.add_argument("--output", action='store', type=str,
        help="Prefix for output files")
    parser.add_argument("--dump", action='store_true',
        help="Dump intermediate files to JSON for troubleshooting")
    args = parser.parse_args()

    # import data
    features = read_pogigwasc_gff(args.features)
    introns = read_realtrons_gff(args.introns)
    # create parent gene features
    genes = create_parent_genes(features)
    # update coordinates of intron-minus gene predictions after introns re-inserted
    intron_parents = update_feature_coords_realtrons(features, introns)
    orphan_introns, conflict_features = check_strand_conflict_orphan_introns(features, introns, intron_parents)
    # segment features that are interrupted by introns
    segment_cds_with_introns(genes, features, introns, intron_parents, orphan_introns, conflict_features)
    # write data in GFF format
    out = []
    for seqid in genes:
        for geneid in genes[seqid]:
            gene_start, gene_end, gene_strand = find_gene_limits(genes, seqid, geneid, features)
            gene_start, gene_end = py2gff_coords(gene_start, gene_end)
            out.append([seqid, 'predicted', 'gene', 
                        gene_start, gene_end, '.', gene_strand, '.', 
                        f'ID={geneid}'])
            for feature_id in genes[seqid][geneid]:
                out.extend(dict2gff(features[seqid][feature_id], feature_id, geneid))
    for seqid in introns:
        for intron_id in introns[seqid]:
            out.extend(dict2gff(introns[seqid][intron_id], intron_id))
    # sort by coordinates
    # negative sort by end coordinates so that enclosing features are listed before their components
    out = sorted(out, key=lambda x: (x[0], int(x[3]), -int(x[4])))
    with open(args.output + ".add_realtrons.gff3", "w") as fh:
        for line in out:
            fh.write("\t".join([str(i) for i in line]))
            fh.write("\n")
    if args.dump:
        with open(args.output + '.dump.genes.json', "w") as fh:
            fh.write(json.dumps(genes, indent=2))
        with open(args.output + '.dump.introns.json', "w") as fh:
            fh.write(json.dumps(introns, indent=2))
        with open(args.output + '.dump.intron_parents.json', "w") as fh:
            fh.write(json.dumps(intron_parents, indent=2))
        with open(args.output + '.dump.orphan_introns.json', "w") as fh:
            fh.write(json.dumps(orphan_introns, indent=2))
        with open(args.output + '.dump.conflict_features.json', "w") as fh:
            fh.write(json.dumps(conflict_features, indent=2))
