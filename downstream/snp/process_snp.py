import os

from downstream.snp.load_1000g import get_1000g


def main():
    snp_dir = "/mnt/stripe/bio/experiments/snp"
    if not os.path.exists(snp_dir):
        os.mkdir(snp_dir)
    get_1000g(snp_dir)


if __name__ == '__main__':
    main()
