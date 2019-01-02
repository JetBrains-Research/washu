Current usage:

$ conda env create --name pipeline-chipseq --file environment.yaml
or
$ conda env update --name pipeline-chipseq --file environment.yaml

$ conda activate pipeline-chipseq
$ snakemake --config work_dir=<work dir> genome=<genome build>