import unittest
import os

from scripts.util import lcs
from scripts.util import find_input

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/input'


class FindInputTest(unittest.TestCase):
    def testLCS(self):
        self.assertTrue(
            lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq', 'Broad_CD14_2_input_ENCFF000CCW.fastq') <
            lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq', 'UW_CD14_input_ENCFF001HUV.fastq'))

    def testFindInput(self):
        self.assertEqual('', find_input(TEST_DATA + '/40_donor6_input.bam'))

        self.assertEqual('40_donor6_input.bam', find_input(TEST_DATA + '/37_donor6_k27ac.bam'))
        self.assertEqual('44_donor7_input.bam', find_input(TEST_DATA + '/41_donor7_k27ac.bam'))

        self.assertEqual('jcl320_wt1_gm_input_ctrl.1919_8.R1_mm10.bam',
                         find_input(TEST_DATA + '/jcl320_ko_gm_h3k27ac_chipd_dna.1919_8.R1_mm10.bam'))
        self.assertEqual('jcl320_wt1_gm_input_ctrl.1919_8.R1_mm10.bam',
                         find_input(TEST_DATA + '/jcl320_wt1_gm_bhlhe40_chipd_dna.1919_8.R1_mm10.bam'))


if __name__ == '__main__':
    unittest.main()
