import unittest

from scripts.macs2_find_input import lcs


class Macs2FindInputTest(unittest.TestCase):
    def test_Input(self):
        self.assertTrue(
            lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq', 'Broad_CD14_2_input_ENCFF000CCW.fastq') <
            lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq', 'UW_CD14_input_ENCFF001HUV.fastq'))


if __name__ == '__main__':
    unittest.main()
