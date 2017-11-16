FROM biolabs/test-data

USER root
# Fix default shell for CONDA source activate command.
RUN ln -snf /bin/bash /bin/sh
RUN apt-get update
RUN apt-get install --yes build-essential libgl1-mesa-dev bc unzip quota
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

RUN conda create -q -n samtools --channel bioconda samtools
RUN conda create -q -n bedtools --channel bioconda bedtools
RUN conda create -q -n r --channel r r
RUN conda create -q -n bowtie --channel bioconda bowtie
RUN conda create -q -n java --channel bioconda fastqc

# Install env py3.5
RUN conda create -n py3.5 python=3.5
# seaborn should be >= 0.8
RUN source activate py3.5 &&\
    conda install --channel bioconda bwa bowtie2 star \
    deeptools sra-tools rseg ucsc-bedgraphtobigwig ucsc-bedclip ucsc-bigwigaverageoverbed && \
    conda install --channel conda-forge matplotlib-venn && \
    conda install pandas numpy scikit-learn pytest pytest-pep8 seaborn && \
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

# To prevent problems with Java interfierence just move execuble to emulate module
RUN mkdir /opt/fastqc && mv /opt/conda/envs/java/bin/fastqc /opt/fastqc


COPY ./test/module.sh /opt/