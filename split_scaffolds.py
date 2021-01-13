#!/usr/bin/env python

import re
import argparse
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, 
    help="Fasta file containing contigs and scaffolds")
parser.add_argument("-o", "--output", type=str, default="test_scaffolds_split",
    help="Output filename prefix")
args = parser.parse_args()


def get_negative_coordinates(pos : list, start : int, end : int):
    """
    get coordinates for the complement to a set of regions

    Regions are defined by a list of pairs of tuples. In addition the start and
    end coordinate for the encompassing region must be provided. Assume that
    coordinates are 0-based end-exclusive (python convention).
    """
    pos = sorted(pos, key = lambda x : x[0]) # sort by start coordinate
    segs = [0]
    for (segstart, segend) in pos:
        segs.extend([segstart, segend])
    segs.append(end)
    return(list(zip(segs[::2], segs[1::2])))


def scaffold2contig_coords(seq_dict):
    """Report non-N regions in input set of scaffolds

    If a sequence has no Ns, the the name is passed through unchanged. Otherwise
    a number is appended to the original scaffold name, to make name of each
    segment derived from splitting up a scaffold on Ns.
    
    Parameters
    ----------
    seq_dict : dict
        Dict of Bio.SeqRecord objects representing scaffolds of interest, keyed by
        scaffold name

    Returns
    -------
    dict
        keyed by new segment name. value is a dict with keys `orig, start, end`
        describing segment position on origial scaffold
    """
    new2old = {}
    for ctg in seq_dict:
        # find stretches of Ns
        nn = re.finditer(r"N+", str(seq_dict[ctg].seq))
        nn_spans = [hit.span() for hit in nn]
        if nn_spans:
            pos_spans = get_negative_coordinates(nn_spans, 0, len(seq_dict[ctg]))
            for i in range(len(pos_spans)):
                newctg = f"{ctg}_{str(i)}"
                new2old[newctg] = { 'orig' : ctg, 
                                    'start' : pos_spans[i][0],
                                    'end' : pos_spans[i][1] }
        else:
            # pass through as is
            new2old[ctg] = { 'orig' : ctg, 'start' : 0, 'end' : len(seq_dict[ctg]) } 
    return(new2old)


def report_split_seqs(seq_dict, new2old):
    """Report sequences of contigs split from scaffolds

    Parameters
    ----------
    seq_dict : dict
    new2old : dict
        Output from scaffold2contig_coords()
    """
    seqlist = []
    for ctg_id in new2old:
        orig = new2old[ctg_id]['orig']
        start = new2old[ctg_id]['start']
        end = new2old[ctg_id]['end']
        ctg = SeqRecord(
            seq_dict[orig].seq[start:end],
            id = ctg_id,
            name = ctg_id,
            description = "")
        seqlist.append(ctg)
    return(seqlist)


if __name__ == "__main__": 
    asm = SeqIO.to_dict(SeqIO.parse(args.input, 'fasta'))
    spls = scaffold2contig_coords(asm)
    asm_new = report_split_seqs(asm, spls)
    SeqIO.write(asm_new, f"{args.output}.fasta", "fasta")
    with open(f"{args.output}.json", "w") as fh:
        fh.write(json.dumps(spls, indent=2))
