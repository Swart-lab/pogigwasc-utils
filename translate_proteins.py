#!/usr/bin/env python

import argparse
import pybedtools

from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


if __name__ == "__main__":
    # command line args
    parser = argparse.ArgumentParser(
        description="""
        Concatenate CDSs and translate to protein sequences

         * Uses karyorelict genetic code (27)
         * Stop codon is removed before translation
         * CDSs belonging to same gene have same Parent in GFF attributes
        """)
    parser.add_argument("-a", "--assembly", action='store', help="""
        Genome assembly, contig names must match GFF file
        """)
    parser.add_argument("-t", "--trans_table", action='store', default=27,
        help="Translation table to use. Default 27 (karyorelict)")
    parser.add_argument("-g", "--gff", action='store', help="""
        Input GFF file from Pogigwasc and/or Pogigwasc + Realtrons pipeline
        """)
    parser.add_argument("-o", "--output", action='store', help="""
        Output file name prefix. Output files will be suffixed `.protein.faa`
        and `.cds.fna`
        """)
    args = parser.parse_args()

    # Read assembly file into memory
    asm = SeqIO.to_dict(SeqIO.parse(args.assembly, 'fasta'))

    # Read GFF file
    gff = pybedtools.BedTool(args.gff)

    # List CDS intervals belonging to each parent gene
    cds_by_parent = defaultdict(list)
    for i in gff:
        if i.fields[2] == 'CDS':
            parent = i.attrs['Parent']
            cds_by_parent[parent].append(i)
    print(f"Gene features: {str(len(cds_by_parent))}")

    # Concatenate CDS intervals and translate
    translated = [] # translations
    problematic = [] # combined CDS lengths not a multiple of 3
    cdss = [] # concatenated CDSs
    for gene in cds_by_parent:
        seq = Seq('')
        # Sort by start interval
        for interval in sorted(cds_by_parent[gene], key=lambda x: x.start):
            seq += asm[interval.chrom][interval.start:interval.end].seq
        if interval.strand == '-':
            seq = seq.reverse_complement()
        if len(seq) % 3 != 0:
            problematic.append(gene)
        else:
            rec = SeqRecord(seq[:-3].translate(table=args.trans_table),
                    id=f'{gene}_trans', name=f'{gene}_trans')
            translated.append(rec)
            # CDS includes stop codon
            cds = SeqRecord(seq, id=f'{gene}', name=f'{gene}')
            cdss.append(cds)
    print(f"Translated gene features: {str(len(translated))}")

    SeqIO.write(translated, args.output+'.protein.faa', 'fasta')
    SeqIO.write(cdss, args.output+'.cds.fna', 'fasta')
