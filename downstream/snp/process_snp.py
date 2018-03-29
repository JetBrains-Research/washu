from downstream.snp.call_snp import make_preprocessed_bams_for_donor
from downstream.snp.snp_data import get_snp_path, get_ucsc_path, get_1000g


donors = ["YD1", "YD2", "YD3", "YD4", "YD5", "YD6", "YD7", "YD8", "YD9", "YD10",
          "YD11", "YD12", "YD14", "YD15", "YD16", "YD17", "YD18", "YD19", "YD20", "YD21",
          "OD1", "OD2", "OD3", "OD4", "OD5", "OD6", "OD7", "OD9", "OD10",
          "OD11", "OD12", "OD13", "OD14", "OD15", "OD16", "OD17", "OD18", "OD19", "OD20"]



def main():
    snp_path = get_snp_path()
    get_1000g(snp_path)
    ucsc_path = get_ucsc_path()

    make_preprocessed_bams_for_donor(snp_path, "OD1")


if __name__ == '__main__':
    main()
