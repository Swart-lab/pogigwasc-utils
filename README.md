Utilities for working with pogigwasc
====================================

Scripts to help prepare data for the pogigwasc gene predictor, and to reformat
or process output files. Requires python 3.6+, Biopython

[pogigwasc](https://github.com/Swart-lab/pogigwasc) is a gene predictor for
eukaryotic genomes with ambiguous stop codons, developed by David Vetter.

TODO:
 * How to split and handle hard masks at beginning or end of contigs?
 * How to deal with empirical introns that are exactly at the boundaries of a CDS?
 * Enhancement: split scaffolds into contigs on soft masks or lower case 'n'
 * Estimate model parameters from training data - standalone script
 * Empirical correction of mispredicted introns from RNAseq data
