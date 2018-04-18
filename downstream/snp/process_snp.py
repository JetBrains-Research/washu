import downstream.snp.call_snp as snp
from downstream.snp.snp_data import get_snp_path, get_ucsc_path, get_1000g, get_dbsnp_path

donors = ["YD1", "YD2", "YD3", "YD4", "YD5", "YD6", "YD7", "YD8", "YD9", "YD10",
          "YD11", "YD12", "YD14", "YD15", "YD16", "YD17", "YD18", "YD19", "YD20", "YD21",
          "OD1", "OD2", "OD3", "OD4", "OD5", "OD6", "OD7", "OD9", "OD10",
          "OD11", "OD12", "OD13", "OD14", "OD15", "OD16", "OD17", "OD18", "OD19", "OD20"]


def main():
    snp_path = get_snp_path()
    get_1000g(snp_path)
    ucsc_path = get_ucsc_path()
    dbsnp_path = get_dbsnp_path()

    paths = []

    for donor in donors:
        paths.append(snp.call_donor_variants(snp_path, ucsc_path, dbsnp_path, donor))

    snp.combine_gvcfs(paths, snp_path, ucsc_path)

    snp.call_snp(snp_path, ucsc_path)

    snp.filter_snp(snp_path, ucsc_path)


if __name__ == '__main__':
    main()
