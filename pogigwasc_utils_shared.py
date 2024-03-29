#!/usr/bin/env python

import re
from collections import defaultdict

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


def gff2py_coords(start:int, end:int):
    """Convert GFF-type 1-based end-inclusive coords to python 0-based end-exclusive"""
    # simple issue but really annoying to keep in head
    # TODO: raise error properly
    if end < start:
        print("invalid GFF coordinates because end {str(end)} before start {str(start)}")
    elif end == start:
        print("zero-length feature")
        return(start, end)
    else:
        return(start-1, end)

def py2gff_coords(start:int, end:int):
    # TODO: zero-length features
    if end < start:
        print("invalid coordinates because end {str(end)} before start (str{start})")
    elif end == start:
        # zero length feature
        return(start, end)
    else:
        return(start+1, end)


def read_pogigwasc_gff(file):
    """
    Import Pogigwasc feature table, convert GFF to python coordinates

    Assumptions:
     - GFF is already sorted by contig and start coordinate when reported by
       Pogigwasc
     - Features do not overlap, stop_codon features are not included in their
       adjacent CDS features
     - phase, score fields are undefined, attributes column not provided by
       Pogigwasc

    If Pogigwasc encounters multiple solutions with equal score for a given
    contig, it will report all of them in 'blocks'. Within each block, the
    features are already sorted by start coordinate. We assume that they are
    all equally good scoring and simply keep the first one, discard the rest.

    Arguments
    ---------
    file : str
        Path to GFF file produced by Pogigwasc

    Returns
    -------
    dict
        key1 - contig ID, key2 - feature ID, value is a dict of feature data
        from GFF file with the following keys: seqid source type old_start
        old_end score strand phase start end. start, end equal old_start,
        old_end; initialized here to update coords subsequently if inserting
        introns. start, end, old_start, old_end coordinates all 0-based
        end-exclusive
    """
    counter = defaultdict( # contig
        lambda: defaultdict( # feature_type
            int)) # running count
    features = defaultdict(dict)
    pgw_gff_keys = ['seqid','source','type','old_start','old_end','score','strand','phase']
    feature_code = { 'gene' : 'g', 'CDS' : 'c', 'stop_codon' : 's' }

    with open(file, "r") as fh:
        # prevent backsliding:
        # pogigwasc reports features in sorted order per contig. if multiple
        # equally scoring solutions, pogigwasc reports multiple 'blocks' of
        # annotations for that contig; use first block, ignore rest
        # TODO: remove duplicates, report conflicting alternatives
        noback = {}
        for line in fh:
            if not re.match(r"\#", line):
                spl = line.rstrip().split("\t")
                # pogigwasc doesn't report attributes column, and leaves score and phase undefined
                spl[3], spl[4] = gff2py_coords(int(spl[3]), int(spl[4]))
                [seqid,source,feature_type,start,end,score,strand,phase] = spl
                if seqid in noback and int(start) < noback[seqid]:
                    # skip if this contig has already been processed up to here
                    pass
                else:
                    counter[seqid][feature_type] += 1
                    if feature_type in feature_code:
                        feature_id = seqid + "." + feature_code[feature_type] + str(counter[seqid][feature_type])
                    else:
                        feature_id = seqid + "." + feature_type + str(counter[seqid][feature_type])
                    # feature_id = "_".join([seqid, feature_type, str(counter[seqid][feature_type])])
                    features[seqid][feature_id] = dict(zip(pgw_gff_keys, spl))
                    noback[seqid] = int(end)
    # initialize new start coords
    for seqid in features:
        for feature_id in features[seqid]:
            features[seqid][feature_id]['start'] = int(features[seqid][feature_id]['old_start'])
            features[seqid][feature_id]['end'] = int(features[seqid][feature_id]['old_end'])
            # Set all CDS phases to zero because Pogigwasc only predicts complete CDSs
            if features[seqid][feature_id]['type'] == 'CDS':
                features[seqid][feature_id]['phase'] = 0

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
                spl[3], spl[4] = gff2py_coords(int( spl[3] ), int( spl[4] ))
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


def flatten_list_of_lists(l : list):
    # given a list containing both lists and non-list elements, flatten to list
    # of elements
    out = []
    for i in l:
        if type(i) is list:
            out.extend(i)
        else:
            out.append(i)
    return(out)


def find_gene_limits(genes, seqid, geneid, features):
    starts = []
    ends = []
    strands = []
    for feature_id in genes[seqid][geneid]:
        if features[seqid][feature_id]['strand'] != '.': # skip features with undefined strand
            strands.append(features[seqid][feature_id]['strand'])
            starts.append(features[seqid][feature_id]['start'])
            ends.append(features[seqid][feature_id]['end'])
    starts = [int(i) for i in flatten_list_of_lists( starts )]
    ends = [int(i) for i in flatten_list_of_lists( ends )]
    strands = set(strands)
    if len(strands) > 1:
        print(f'strand conflict within component features of gene {geneid}')
    else:
        if min(starts) > min(ends) or max(ends) < max(starts):
            print(f"start before end in component features of gene {geneid}")
        else:
            return(min(starts), max(ends), list(strands)[0])


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


