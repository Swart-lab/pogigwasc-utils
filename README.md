Utilities for working with pogigwasc
====================================

Scripts to help prepare data for the pogigwasc gene predictor, and to reformat
or process output files, and to work with output from 
[Intronarrator](https://github.com/Swart-lab/Intronarrator).

Requires python 3.6+, Biopython, pybedtools

[pogigwasc](https://github.com/Swart-lab/pogigwasc) is a gene predictor for
eukaryotic genomes with ambiguous stop codons, developed by David Vetter.

Scripts
-------

 * `split_scaffolds.py` - Split scaffolds on Ns into contigs, and record
   original coordinates in a JSON file (used later to recombine the contigs).
   This is because pogigwasc does accepts only contigs (no Ns) as input
 * `recombine_contigs.py` - Recombine annotations performed on split contigs to
   original scaffold coordinates. Requires JSON output from `split_scaffolds.py`
 * `add_attributes_pogigwasc.py` - Add parent gene features and ID, Parent
   attributes to pogigwasc predictions and extend CDSs to include adjacent stop
   codons.
 * `add_realtrons_pogigwasc_intronless.py` - Combine pogigwasc gene predictions
   and "realtrons" (empirically predicted introns from RNAseq mapping) predicted
   by Intronarrator. Intronarrator must first be used to identify realtrons and
   artificially remove them from the assembly. Gene prediction then run with
   pogigwasc on the intronless assembly in intronless mode.
 * `pogigwasc_utils_shared.py` - Functions shared by more than one script (not
   used as a standalone script)
 * `testing.py` - Test module for these scripts


TODO:
 * How to split and handle hard masks at beginning or end of contigs?
 * How to deal with empirical introns that are exactly at the boundaries of a CDS?
 * Enhancement: split scaffolds into contigs on soft masks or lower case 'n'
 * Estimate model parameters from training data - standalone script
 * Empirical correction of mispredicted introns from RNAseq data
