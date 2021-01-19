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


def renumber_features(new2old, gff_file):
    out = []
    with open(gff_file, "r") as fh:
        for line in fh:
            if re.match(r"\#", line):
                out.append([line.rstrip()])
            else:
                ll = line.rstrip().split("\t")
                # seqid source type start end score strand phase attributes
                try:
                    scaff = new2old[ll[0]]["orig"]
                    # no +/- 1 required
                    # output is in 1-based end-inclusive coordinates
                    scaff_start = int(ll[3]) + int(new2old[ll[0]]["start"])
                    scaff_end = int(ll[4]) + int(new2old[ll[0]]["start"])
                    out.append([scaff, ll[1], ll[2], scaff_start, scaff_end] + ll[5:])
                except KeyError:
                    print(f"contig {ll[0]} not found in JSON file, please check")
    return(out)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
        Recombine split contigs into scaffolds, adding Ns to fill gaps
        """)
    parser.add_argument( # compulsory
        "-j", "--json", action='store',
        help="JSON file containing mapping of coordinates from split contigs to original scaffolds")
    parser.add_argument( # optional
        "-c", "--contigs", action='store',
        help="Fasta file containing split contigs")
    parser.add_argument( # optional
        "-f", "--features", action='store',
        help="Feature table in GFF format for contigs, with coordinates to be updated to recombined scaffolds")
    parser.add_argument(
        "-o", "--output", action="store",
        help="Output filename prefix to write recombined scaffolds and/or modified feature table")

    args = parser.parse_args()

    with open(args.json, "r") as fh:
        new2old = json.load(fh)
    old2new = rekey_regions_by_original_contig(new2old) 
    # Report updated contigs
    if args.contigs:
        asm = SeqIO.to_dict(SeqIO.parse(args.contigs, 'fasta'))
        # assume that scaffolds do not start or end with Ns! 
        scaff_asm = add_gaps(old2new, asm)
        SeqIO.write([scaff_asm[scaff] for scaff in scaff_asm], args.output+'.scaffolds.fasta', 'fasta')
    if args.features:
        gff_out = renumber_features(new2old, args.features)
        with open(args.output+'.features.gff3', "w") as fh:
            for ll in gff_out:
                fh.write("\t".join([str(i) for i in ll]) + "\n")
