#!/usr/bin/env python

import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, 
    help="Fasta file containing contigs and scaffolds")
parser.add_argument("-o", "--output", type=str, default="test_scaffolds_split",
    help="Output filename prefix")
parser.add_argument("-w", "--width", type=int, default=80,
    help="Fasta file output width")
args = parser.parse_args()

def wrap_line_buffer(line, fh, width=80):
    """Wrap a line to desired width, write full width lines, and return remainder
    that does not fill complete line"""
    if len(line) > width:
        for rep in range(int(len(line)/width)): # floor line length / width
            fh.write(line[rep*width : (rep+1)*width] + "\n")
        return(line[int(len(line)/width)*width:])
    else:
        return(line)

tbl = []
# list of lists, each with three fields
# new contig id / old contig id / coordinate offset of new vs old contig

fh_in = open(args.input, "r")
fh_out = open(args.output + ".fasta", "w")

ctg = None
newctg = None
linebuff = ""
coord = 0
for line in fh_in:
    line = line.rstrip()
    if re.match(">", line):
        # Dump data from old contig
        if ctg:
            fh_out.write(linebuff + "\n")
            linebuff = "" # Clear buffer
            coord = 0
        # Start new contig
        ctg = line[1:] # strip starting > character
        splitcount = 0
        fh_out.write(">" + ctg + "\n")
    else:
        ns = re.finditer(r"N+", line)
        ncoords = [i.span() for i in ns]
        if len(ncoords) > 0:
            # Remove Ns and record coordinates
            # Add to line buffer
            start = 0
            for (nstart, nend) in ncoords:
                linebuff += line[start:nstart]
                linebuff = wrap_line_buffer(linebuff, fh_out, args.width)
                if len(linebuff) > 0:
                    fh_out.write(linebuff + "\n")
                linebuff = ""
                if nend < len(line):
                    # Start a new record
                    splitcount += 1
                    newctg = ctg + "_" + str(splitcount)
                    fh_out.write(">" + newctg + "\n")
                    newcoord = coord + nend
                    tbl.append([newctg, ctg, newcoord])
                # Update start coordinate for next non-N segment
                start = nend
            if start < len(line):
                # catch trailing sequence that is not N
                linebuff += line[start:]
                linebuff = wrap_line_buffer(linebuff, fh_out, args.width)
            coord += len(line)
        else:
            # Add to existing line
            linebuff += line
            coord += len(line)
            linebuff = wrap_line_buffer(linebuff, fh_out, args.width)
# Print trailing buffer sequence
if len(linebuff) > 0:
    fh_out.write(linebuff + "\n")
# Close file handles
fh_in.close()
fh_out.close()

# Write table of coordinate offsets for new contigs
with open(args.output + "table.tsv", "w") as fh:
    for l in tbl:
        fh.write("\t".join([str(i) for i in l]) + "\n")
