#!/usr/bin/env python3

import re
from collections import defaultdict
import argparse

GFF_KEYS = ['seqid','source','type','start','end','score','strand','phase','attributes']


def length2phase(l: int, p=0):
    """Convert nucleotide sequence length to phase
    
    Given a base in a CDS that is l bases from the start, what is
    the phase of that base? 
    
    From the GFF3 specification:
    > The phase is one of the integers 0, 1, or 2, indicating the number of
    > bases forward from the start of the current CDS feature the next codon
    > begins.
    
    This is different from frame!
    
    Parameters
    ----------
    l : int
        number of bases from start (or frame-0 base) to the base for which the
        phase should be calculated. This includes the start base but excludes
        the base of interest.
    p : int
        phase of the first base; default 0, meaning that the preceding segment
        is from the start of the CDS.
    """
    f2p = { 0:0, 1:2, 2:1 } # convert frame to phase
    # convert initial phase to frame
    try:
        f = f2p[p]
    except KeyError:
        print(f"invalid phase {str(p)} specified, should be 0,1,2")
    
    frame = (l + f) % 3
    return(f2p[frame])


def read_pogigwasc_gff(file):
    # import minus_introns feature table
    counter = defaultdict( # contig
        lambda: defaultdict( # feature_type
            int)) # running count
    features = defaultdict(dict)
    pgw_gff_keys = ['seqid','source','type','old_start','old_end','score','strand','phase']

    with open(file, "r") as fh:
        for line in fh:
            if not re.match("\#", line):
                spl = line.rstrip().split("\t")
                # pogigwasc doesn't report attributes column, and leaves score and phase undefined
                [seqid,source,feature_type,start,end,score,strand,phase] = spl
                counter[seqid][feature_type] += 1
                feature_id = "_".join([seqid, feature_type, str(counter[seqid][feature_type])])
                features[seqid][feature_id] = dict(zip(pgw_gff_keys,spl))
    # initialize new start coords
    for seqid in features:
        for feature_id in features[seqid]:
            features[seqid][feature_id]['start'] = int(features[seqid][feature_id]['old_start'])
            features[seqid][feature_id]['end'] = int(features[seqid][feature_id]['old_end'])

    return(features)


def read_realtrons_gff(file):
    # import minus_introns feature table
    counter = defaultdict( # contig
        lambda: defaultdict( # feature_type
            int)) # running count
    # import introns table
    introns = defaultdict(dict)
    with open(file, "r") as fh:
        for line in fh:
            if not re.match("\#", line):
                spl = line.rstrip().split("\t")
                [seqid,source,feature_type,start,end,score,strand,phase,attributes] = spl
                attrs = attributes.rstrip(';').split(';')
                attrdict = {a.split('=')[0] : a.split('=')[1] for a in attrs}
                counter[seqid][feature_type] += 1
                if 'ID' in attrdict:
                    intron_id = attrdict['ID']
                else:
                    intron_id = "_".join([seqid, feature_type, str(counter[seqid][feature_type])])
                introns[seqid][intron_id] = dict(zip(GFF_KEYS, spl))
                introns[seqid][intron_id]['start'] = int(introns[seqid][intron_id]['start'])
                introns[seqid][intron_id]['end'] = int(introns[seqid][intron_id]['end'])
    return(introns)


def update_feature_coords_realtrons(features, introns):
    # Adjust coordinates of features with introns added
    # Record features that contain introns
    # features is changed in-place

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
                intron_len = introns[seqid][intron_id]['end'] - introns[seqid][intron_id]['start'] + 1 # plus 1 because end-inclusive

                for feature_id in feature_ids_sorted:
                    if features[seqid][feature_id]['start'] > introns[seqid][intron_id]['start']:
                        # feature outside of intron: push whole feature forward
                        features[seqid][feature_id]['start'] += intron_len
                        features[seqid][feature_id]['end'] += intron_len
                    elif features[seqid][feature_id]['end'] >= introns[seqid][intron_id]['start']:
                        # start of intron within un-adjusted coordinates of feature
                        # keep start coordinate, push end coordinate forward
                        features[seqid][feature_id]['end'] += intron_len
                        intron_parents[seqid][feature_id].append(intron_id)

        except KeyError:
            print(f"seqid {seqid} in intron gff not found in feature gff")
            pass
    return(intron_parents)


def check_strand_conflict_orphan_introns(features, introns, intron_parents):
    # Report adjusted coordinates

    orphan_introns = defaultdict(list)
    conflict_features = defaultdict( # seqid
        lambda: defaultdict( # feature id
            list))

    for seqid in features:
        # sort introns ascending by start coordinate
        introns_list = [ (i, introns[seqid][i]) for i in introns[seqid]]
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
    # assumes that every CDS is adjacent to a stop codon,
    # and that a gene is composed of no more than one CDS+stop codon

    # Run this after update_feature_coords_realtrons()
    # modifies features object in-place

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
                if int(features[seqid][last_feature]['end']) + 1 == int(features[seqid][feature_id]['start']):
                    # either: CDS followed by stop codon and strands both positive
                    # or: stop codon followed by CDS, and strands both negative
                    if (features[seqid][last_feature]['type'] == 'CDS' \
                            and features[seqid][feature_id]['type'] == 'stop_codon' \
                            and features[seqid][last_feature]['strand'] == '+' \
                            and features[seqid][feature_id]['strand'] == '+') \
                        or (features[seqid][last_feature]['type'] == 'stop_codon' \
                            and features[seqid][feature_id]['type'] == 'CDS' \
                            and features[seqid][last_feature]['strand'] == '-' \
                            and features[seqid][feature_id]['strand'] == '-'):
                        # update gene counter
                        geneid_counter[seqid] += 1
                        geneid = "_".join([seqid, 'gene', str(geneid_counter[seqid])])
                        # add feature ID and Parent gene id to attributes field
                        genes[seqid][geneid] = { feature_id : features[seqid][feature_id],
                                                 last_feature: features[seqid][last_feature]}
                    else:
                        print("adjacent features but apparently unrelated")
                # Update previous feature
                last_feature = feature_id
            else:
                # first entry
                last_feature = feature_id
    return(genes)


def segment_cds_with_introns(genes, introns, intron_parents, conflict_features):
    out = []
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
                        genes[seqid][geneid][feature_id]['Note'] = 'strand_conflict'
                        # report overlapping introns as 'intron_orphan' features
                        for conflict_intron in conflict_features[seqid][feature_id]:
                            introns[seqid][conflict_intron]['Note'] = 'orphan'

                    else:
                        # check if CDS feature (or start_codon or stop_codon) and split if necessary
                        if genes[seqid][geneid][feature_id]['type'] in ['CDS','start_codon','stop_codon']:
                            # get coordinates of the split CDS feature segments,
                            # interrupted by introns
                            seg_coords = [genes[seqid][geneid][feature_id]['start']]
                            for intron_id in intron_parents[seqid][feature_id]:
                                # the introns should be present in coordinate sort order
                                # +/- 1 because GFF coordinates are end-inclusive
                                seg_coords.append(introns[seqid][intron_id]['start']-1)
                                seg_coords.append(introns[seqid][intron_id]['end']+1)
                                # add gene as parent to this intron
                                introns[seqid][intron_id]['Parent'] = geneid
                            seg_coords.append(features[seqid][feature_id]['end'])

                            # start and end coords as a list for each segment
                            starts = seg_coords[::2]
                            ends = seg_coords[1::2]
                            segs = []
                            # forward strand
                            if genes[seqid][geneid][feature_id]['strand'] == '+':
                                segs = list(zip(starts,ends))
                            # reverse strand: work backwards
                            elif genes[seqid][geneid][feature_id]['strand'] == '-':
                                segs = reversed(list(zip(starts,ends)))
                            # get list of phases
                            # initialize phase, each intronless CDS starts from 0
                            phase = 0
                            phases = []
                            for (i,j) in segs:
                                seglen = j - i + 1
                                phases.append(phase)
                                #update phase for next segment
                                phase = length2phase(seglen, phase)
                            # reverse because segments are iterated in reverse order
                            if genes[seqid][geneid][feature_id]['strand'] == '-':
                                phases = list(reversed(phases))
                            # replace feature 'start', 'end', and 'phase' fields with lists
                            genes[seqid][geneid][feature_id]['start'] = starts
                            genes[seqid][geneid][feature_id]['end'] = ends
                            genes[seqid][geneid][feature_id]['phase'] = phases
                        else:
                            print(f"Non CDS feature {feature_id} overlaps with empirical intron")
                            # this should never happen because update_feature_coords_realtrons() 
                            # is called with pogigwasc feature table which has only CDS and stop_codons
                else:
                    # intronless feature, report as-is, change phase from '.' to 0 if CDS
                    if genes[seqid][geneid][feature_id]['type'] == 'CDS':
                        genes[seqid][geneid][feature_id]['phase'] = 0
        # add orphan introns
        for orphan_id in orphan_introns[seqid]:
            introns[seqid][orphan_id]['Note'] = 'orphan'


def dict2gff(feature: dict, feature_id, parent_id=None):
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
            line = [feature['seqid'], feature['source'], feature['type'],
                    feature['start'][i], 
                    feature['end'][i],
                    feature['score'], feature['strand'],
                    feature['phase'][i],
                    attributes_string]
            out.append(line)

    else:
        line = [feature[col] for col in GFF_KEYS[0:8]]
        line.append(attributes_string)
        out.append(line)
    
    return(out)


def find_gene_limits(gene, geneid):
    starts = []
    ends = []
    strands = []
    for feature_id in gene:
        if gene[feature_id]['strand'] != '.': # skip features with undefined strand
            strands.append(gene[feature_id]['strand'])
        if type(gene[feature_id]['start']) is list:
            starts.extend(gene[feature_id]['start'])
            ends.extend(gene[feature_id]['end'])
        elif type(gene[feature_id]['start']) is int or str:
            starts.append(gene[feature_id]['start'])
            ends.append(gene[feature_id]['start'])
    starts = [int(i) for i in starts]
    ends = [int(i) for i in ends]
    strands = set(strands)
    if len(strands) > 1:
        print(f'strand conflict within component features of gene {geneid}')
    else:
        if min(starts) > min(ends) or max(ends) < max(starts):
            print(f"start before end in component features of gene {geneid}")
        else:
            return(min(starts), max(ends), list(strands)[0])


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--introns", action='store', type=str, 
        help="Realtrons GFF3 file")
    parser.add_argument("--features", action='store', type=str, 
        help="Pogigwasc intronless-mode gene predictions, GFF3 file")
    parser.add_argument("--output", action='store', type=str,
        help="Prefix for output files")
    args = parser.parse_args()

    # import data
    features = read_pogigwasc_gff(args.features)
    introns = read_realtrons_gff(args.introns)
    # update coordinates of intron-minus gene predictions after introns re-inserted
    intron_parents = update_feature_coords_realtrons(features, introns)
    orphan_introns, conflict_features = check_strand_conflict_orphan_introns(features, introns, intron_parents)
    # create parent gene features
    genes = create_parent_genes(features)
    # segment features that are interrupted by introns
    segment_cds_with_introns(genes, introns, intron_parents, conflict_features)
    # write data in GFF format
    out = []
    for seqid in genes:
        for geneid in genes[seqid]:
            gene_start, gene_end, gene_strand = find_gene_limits(genes[seqid][geneid], geneid)
            out.append([seqid, 'prediction', 'gene', 
                        gene_start, gene_end, '.', gene_strand, '.', 
                        f'ID={geneid}'])
            for feature_id in genes[seqid][geneid]:
                out.extend(dict2gff(genes[seqid][geneid][feature_id], feature_id, geneid))
    for seqid in introns:
        for intron_id in introns[seqid]:
            out.extend(dict2gff(introns[seqid][intron_id], intron_id))
    # sort by coordinates
    # negative sort by end coordinates so that enclosing features are listed before their components
    out = sorted(out, key=lambda x: (int(x[3]), -int(x[4])))
    with open(args.output + ".add_realtrons.gff3", "w") as fh:
        for line in out:
            fh.write("\t".join([str(i) for i in line]))
            fh.write("\n")
