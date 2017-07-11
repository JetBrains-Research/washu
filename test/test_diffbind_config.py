import os
import unittest

from bed.bedtrace import run

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + '/testdata/diffbind_config'
DIFFBIND_CONFIG_SH = os.path.dirname(os.path.abspath(__file__)) + '/../analysis/diffbind_config.sh'


class DiffBindConfigTest(unittest.TestCase):

    def test_diffbind_config(self):
        config = run(
            [['bash', DIFFBIND_CONFIG_SH, TEST_DATA + '/k4me3_bams', TEST_DATA + '/k4me3_bams_macs2_broad_0.1']]
        )[0].decode('utf-8').replace(TEST_DATA + '/', '')
        print(config)
        self.assertEqual("""SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
OD1,CD14,Age,O,1,k4me3_bams/OD1_k4me3_hg19.bam,O_input,OD_input.bam,k4me3_bams_macs2_broad_0.1/OD1_k4me3_hg19_broad_0.1_peaks.xls,macs
""", config)


if __name__ == '__main__':
    unittest.main()
