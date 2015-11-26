from __future__ import print_function
import os
import unittest
import operator
import tempfile
import shutil
import logging

import rabifier
import rabifier.rabmyfire


logging.basicConfig(level=logging.WARNING)

SEQUENCES = os.path.join(os.path.dirname(__file__), 'data/sequences.fa')


class Phase1Test(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.phase1 = rabifier.rabmyfire.Phase1(SEQUENCES, self.tmpdir, cpu=4)

    def tearDown(self):
        self.phase1 = None
        shutil.rmtree(self.tmpdir)

    def test_gdomain(self):
        protein2gdomain = self.phase1.find_gdomain()
        self.assertEqual(protein2gdomain['rab8'][0], '8-118')
        self.assertEqual(protein2gdomain['rab11'][0], '6-188')
        self.assertNotIn('myoglobin', protein2gdomain)

    def test_best_hit(self):
        best_hit = self.phase1.find_best_hit()
        self.assertAlmostEqual(best_hit['rab8']['evalue_bh_rabs'], 7.3e-73)
        self.assertAlmostEqual(best_hit['rab8']['evalue_bh_non_rabs'], 1.2e-28)
        self.assertAlmostEqual(best_hit['rab11']['evalue_bh_rabs'], 1.1e-149)
        self.assertAlmostEqual(best_hit['rab11']['evalue_bh_non_rabs'], 3.9e-43)
        self.assertIsNone(best_hit['myoglobin']['evalue_bh_rabs'])
        self.assertIsNone(best_hit['myoglobin']['evalue_bh_non_rabs'])

    def test_rabf_motif_scanner(self):
        rabf = self.phase1.find_rabf_motifs()
        self.assertListEqual(list(map(operator.itemgetter(0), rabf['rab8'])), [1, 2, 3, 4, 5])
        self.assertListEqual(list(map(operator.itemgetter(0), rabf['rab11'])), [1, 3, 4, 5])
        self.assertListEqual(rabf['myoglobin'], [])


class Phase2Test(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_subfamily_classification(self):
        phase2 = rabifier.rabmyfire.Phase2(SEQUENCES, self.tmpdir)
        phase2()
        self.assertAlmostEqual(phase2.score['rab8']['rab15'], 0.130272143641)
        self.assertAlmostEqual(phase2.score['rab11']['rab11'], 0.96194485463)
        self.assertEqual(phase2.score['myoglobin']['rabX'], 0.5)


if __name__ == '__main__':
    unittest.main()
