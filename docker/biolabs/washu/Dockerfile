FROM biolabs/test-data

USER root
# Fix default shell for CONDA source activate command.
RUN ln -snf /bin/bash /bin/sh
RUN apt-get update

# Essential development tools and headers
RUN apt-get install --yes build-essential libgl1-mesa-dev bc unzip quota libssl-dev \
libboost-dev curl libssh2-1-dev libcurl4-openssl-dev

# GNU AWK requered for proper scripts work.
RUN apt-get install --yes gawk

# Hack to enable MACS2 in another conda environment
RUN conda create -n macs2 --channel bioconda macs2
RUN ln -sf /opt/conda/envs/macs2/bin/macs2 /usr/local/bin/macs2

# SICER
RUN pip install scipy
RUN cd /tmp && wget http://home.gwu.edu/~wpeng/SICER_V1.1.tgz && tar xvf SICER_V1.1.tgz && mv SICER_V1.1 /opt
# Please refer to README for installation instructions, modify scripts, i.e.
RUN sed -i 's#/home/data/SICER1.1#/opt/SICER_V1.1#g' /opt/SICER_V1.1/SICER/SICER.sh
# SICER is python2 library, force it!
RUN sed -i 's#python#python2#g' /opt/SICER_V1.1/SICER/SICER.sh
RUN chmod a+x /opt/SICER_V1.1/SICER/SICER.sh
ENV PATH $PATH:/opt/SICER_V1.1/SICER

# Install packages to be loaded using module in dedicated envs
RUN conda create -y -q -n samtools --channel bioconda samtools
RUN conda create -y -q -n bedtools --channel bioconda bedtools
RUN conda create -y -q -n r --channel r r
RUN conda create -y -q -n bowtie --channel bioconda bowtie
RUN conda create -y -q -n bowtie2 --channel bioconda bowtie2
RUN conda create -y -q -n java --channel bioconda fastqc
RUN conda create -y -q -n sratoolkit --channel bioconda sra-tools

# Deeptools installs samtools and lots of dependencies, isolate it
RUN conda create -y -q -n deeptools --channel bioconda deeptools
RUN ln -s /opt/conda/envs/deeptools/bin/deeptools /usr/local/bin/deeptools
RUN ln -s /opt/conda/envs/deeptools/bin/bamCoverage /usr/local/bin/bamCoverage

# Install env py3.5
RUN conda create -y -n py3.5 python=3.5
# seaborn should be >= 0.8
RUN source activate py3.5 && \
    conda install -y --channel bioconda star rseg \
    ucsc-bedgraphtobigwig ucsc-bedclip ucsc-bigwigaverageoverbed && \
    conda install -y --channel conda-forge matplotlib-venn && \
    conda install -y pandas numpy scikit-learn pytest pytest-pep8 seaborn lxml && \
    pip install multiqc teamcity-messages

# Workaround for TeamCity CI, temp folders are created with root permissions, unacessible for USER
# RUN groupadd -r washu && useradd -ms /bin/bash -g washu user
# WORKDIR /home/user
# USER user

# Download Picard tools
RUN cd ~ && wget -q https://github.com/broadinstitute/picard/releases/download/2.10.7/picard.jar

# Download ZINBRA
RUN cd ~ && wget -q https://github.com/JetBrains-Research/zinbra/releases/download/v0.4.0/zinbra-0.4.0.jar -O zinbra.jar
# Alternative CI link
# https://teamcity.jetbrains.com/repository/download/Epigenome_Zinbra/lastPinned/zinbra-0.4.0.jar?guest=1

# Install phantom peak qual tools
RUN cd ~ && wget -q https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz
RUN cd ~ && tar xvf ccQualityControl.v.1.1.tar.gz

# Install spp R package for phantom peak qual tools
RUN /opt/conda/envs/r/bin/R -e "source('https://bioconductor.org/biocLite.R');biocLite(c('caTools', 'Rsamtools'));\
install.packages('devtools', repos='http://cran.us.r-project.org');\require(devtools);\
devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)"

# To prevent problems with Java interfierence just move execuble to emulate module
RUN mkdir /opt/fastqc && mv /opt/conda/envs/java/bin/fastqc /opt/fastqc


COPY ./module.sh /opt/