#!/usr/bin/env python

import re
import unittest
from sys import argv
import os.path

# Prefix to folder where testing.py is located
# subfolder test/ containing test files should be in the same folder
prefix = os.path.abspath(os.path.dirname(argv[0]))

from add_realtrons_pogigwasc_intronless import *

# Run test: `python -v testing.py`
# with -v switch for verbose mode

class Test_add_realtrons_pogigwasc_intronless(unittest.TestCase):
    
    # utility functions
    def test_length2phase(self):
        self.assertEqual( length2phase(0,0), 0)
        self.assertEqual( length2phase(0,2), 2)
        self.assertEqual( length2phase(0,1), 1)
        # multiples of 3 do not change phase
        self.assertEqual( length2phase(3,0), 0)
        self.assertEqual( length2phase(9,1), 1)
        self.assertEqual( length2phase(3,2), 2)
        # frame of +1 shifts phase by 2
        self.assertEqual( length2phase(4,0), 2)
        self.assertEqual( length2phase(4,1), 0)
        self.assertEqual( length2phase(4,2), 1)

    def test_coordinate_conversion(self):
        self.assertEqual( gff2py_coords(15, 100),
            (14, 100))
        self.assertEqual( py2gff_coords(14, 100),
            (15, 100))

    def test_flatten_list_of_lists(self):
        lol = [[12, 35, 66], 80, [93, 120, 150], 200, 300]
        flat = [12, 35, 66, 80, 93, 120, 150, 200, 300]
        self.assertEqual( flatten_list_of_lists(lol), flat)

    def test_create_parent_genes(self):
        # read GFF and convert coords to python type
        # original: (47,70), (71,73)
        features = read_pogigwasc_gff(os.path.join(prefix,'test/test_pogigwasc_intronless.gff3'))
        self.assertEqual( features['contig_1']['contig_1_CDS_1']['start'], 46 )
        self.assertEqual( features['contig_1']['contig_1_CDS_1']['end'], 70 )
        self.assertEqual( features['contig_1']['contig_1_stop_codon_1']['start'], 70 )
        self.assertEqual( features['contig_1']['contig_1_stop_codon_1']['end'], 73 )
        # create parent gene features and extend CDS to encompass stop codon
        genes = create_parent_genes(features)
        self.assertTrue( 'contig_1_CDS_1' in genes['contig_1']['contig_1_gene_1'])
        self.assertTrue( 'contig_1_stop_codon_1' in genes['contig_1']['contig_1_gene_1'])
        self.assertTrue( 'contig_1_CDS_2' in genes['contig_1']['contig_1_gene_2'])
        self.assertTrue( 'contig_1_stop_codon_2' in genes['contig_1']['contig_1_gene_2'])

    def test_add_introns_to_features(self):
        # read features
        features = read_pogigwasc_gff(os.path.join(prefix,'test/test_pogigwasc_intronless.gff3'))
        genes = create_parent_genes(features)
        # read introns 
        introns = read_realtrons_gff(os.path.join(prefix,'test/test_realtrons.gff3'))

        # identify parent features for introns and extend features
        intron_parents = update_feature_coords_realtrons(features, introns)
        self.assertEqual( intron_parents['contig_1']['contig_1_CDS_1'],
            ['intron1', 'intron1a', 'intron1b'])
        self.assertEqual( intron_parents['contig_1']['contig_1_stop_codon_1'],
            ['intron1b'])
        self.assertEqual ( features['contig_1']['contig_1_CDS_1']['start'], 46)
        self.assertEqual ( features['contig_1']['contig_1_CDS_1']['end'], 120)

        # identify orphan introns and conflicts
        orphan_introns, conflict_features = check_strand_conflict_orphan_introns(
            features, introns, intron_parents)
        self.assertTrue ('intron_orphan' in orphan_introns['contig_1'])
        # TODO: make test for features/introns with conflicting strands

        # split CDS features that overlap with introns
        segment_cds_with_introns(genes, features, introns,
            intron_parents, orphan_introns, conflict_features)
        self.assertEqual( features['contig_1']['contig_1_CDS_1']['start'],
            [46, 67, 83, 118])
        self.assertEqual( features['contig_1']['contig_1_CDS_1']['end'],
            [46, 73, 102, 120])
        self.assertEqual( features['contig_1']['contig_1_stop_codon_1']['start'],
            [101, 118])
        self.assertEqual( features['contig_1']['contig_1_stop_codon_1']['end'],
            [102, 120])
        self.assertEqual( introns['contig_1']['intron1']['Parent'],
            'contig_1_gene_1')
        self.assertEqual( introns['contig_1']['intron1a']['Parent'],
            'contig_1_gene_1')
        self.assertEqual( introns['contig_1']['intron1b']['Parent'],
            'contig_1_gene_1')

    def test_find_gene_limits(self):
        # read features
        features = read_pogigwasc_gff(os.path.join(prefix,'test/test_pogigwasc_intronless.gff3'))
        genes = create_parent_genes(features)
        # read introns 
        introns = read_realtrons_gff(os.path.join(prefix,'test/test_realtrons.gff3'))
        # identify parent features for introns and extend features
        intron_parents = update_feature_coords_realtrons(features, introns)
        # identify orphan introns and conflicts
        orphan_introns, conflict_features = check_strand_conflict_orphan_introns(
            features, introns, intron_parents)
        # split CDS features that overlap with introns
        segment_cds_with_introns(genes, features, introns,
            intron_parents, orphan_introns, conflict_features)
        self.assertEqual(
            find_gene_limits(genes, 'contig_1', 'contig_1_gene_1', features),
            (46, 120, '+'))


if __name__ == "__main__":
    unittest.main()
