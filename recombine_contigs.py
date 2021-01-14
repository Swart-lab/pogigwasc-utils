#!/usr/bin/env python3

import re
import json
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def complement_coordinates_inside(pos : list):
    # look for internal complementary regions
    pos = sorted(pos, key = lambda x : x[0]) # sort by start coordinate
    segs = []
    for i in pos:
        # assume that each tuple contains start,end coordinates only
        segs.extend([i[0],i[1]])
    segs = segs[1:len(segs)-1] # remove first and last values
    return(list(zip(segs[::2], segs[1::2])))


def rekey_regions_by_original_contig(new2old : dict):
    old2new = defaultdict(list)
    for ctg in new2old:
        scaff = new2old[ctg]['orig']
        old2new[scaff].append(
            (ctg, new2old[ctg]['start'], new2old[ctg]['end'], 'contig')
            )
    return(old2new)


def add_gaps(old2new, asm):
    # change old2new in place
    # get coordinates of N regions and create dummy sequences for them
    newasm = {}
    for scaff in old2new:
        ctg_coords = [(i[1],i[2]) for i in old2new[scaff]]
        ctg_coords = sorted(ctg_coords, key=lambda x: x[0]) # sort by start coord
        gaps = complement_coordinates_inside(ctg_coords) # get gap coords
        gap_names = ['gap_' + str(i) for i in range(len(gaps))] 
        # assume that real contigs are not named 'gap_' !
        gaps = [(gap_names[i], gaps[i][0], gaps[i][1], 'gap') for i in range(len(gaps))]
        old2new[scaff].extend(gaps) # add gaps
        old2new[scaff] = sorted(old2new[scaff], key=lambda x: x[1]) # sort by start
        # create gapped sequences
        scaffseq = "" # string
        for (ctg, start, end, seqtype) in old2new[scaff]:
            if seqtype == 'contig':
                # get ctg sequence
                scaffseq = scaffseq + str(asm[ctg].seq)
            elif seqtype == 'gap':
                # create gap sequence
                gapseq = 'N' * (end - start)
                scaffseq = scaffseq + gapseq
        newasm[scaff] = SeqRecord(
            Seq(scaffseq), name=scaff, id=scaff, description="")
    return(newasm)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
        Recombine split contigs into scaffolds, adding Ns to fill gaps
        """)
    parser.add_argument(
        "-s", "--scaffolds", action='store',
        help="Fasta file containing split contigs")
    parser.add_argument(
        "-o", "--output", action="store",
        help="Path to output file to write recombined scaffolds")
    parser.add_argument(
        "-j", "--json", action='store',
        help="JSON file containing mapping of coordinates from split contigs to original scaffolds")

    args = parser.parse_args()

    asm = SeqIO.to_dict(SeqIO.parse(args.scaffolds, 'fasta'))
    with open(args.json, "r") as fh:
        spls = json.load(fh)
    # assume that scaffolds do not start or end with Ns! 
    old2new = rekey_regions_by_original_contig(spls)
    scaff_asm = add_gaps(old2new, asm)
    SeqIO.write([scaff_asm[scaff] for scaff in scaff_asm], args.output, 'fasta')
