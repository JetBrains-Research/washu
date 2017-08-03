FROM continuumio/miniconda

USER root
# Fix default shell for CONDA source activate command.
RUN ln -snf /bin/bash /bin/sh
RUN apt-get update
RUN apt-get install --yes build-essential libgl1-mesa-dev bc

# Hack to enable MACS2 in another conda environment
RUN conda install --channel bioconda macs2
RUN ln -sf $(which macs2) /usr/local/bin/macs2

# SICER
RUN cd /tmp && wget http://home.gwu.edu/~wpeng/SICER_V1.1.tgz && tar xvf SICER_V1.1.tgz && mv SICER_V1.1 /opt
# Please refer to README for installation instructions, modify scripts, i.e.
RUN sed -i 's#/home/data/SICER1.1#/opt/SICER_V1.1#g' /opt/SICER_V1.1/SICER/SICER.sh
# SICER is python2 library, force it!
RUN sed -i 's#python#python2#g' /opt/SICER_V1.1/SICER/SICER.sh
RUN chmod a+x /opt/SICER_V1.1/SICER/SICER.sh
ENV PATH $PATH:/opt/SICER_V1.1/SICER

# Install env py3.5
RUN conda create -n py3.5 python=3.5
RUN source activate py3.5 &&\
    conda install --channel bioconda fastqc bwa bowtie bowtie2 star samtools bedtools \
    deeptools sra-tools rseg ucsc-bedgraphtobigwig ucsc-bedclip &&\
    conda install --channel conda-forge matplotlib-venn && conda install --channel r r && conda install pandas numpy &&\
    pip install multiqc teamcity-messages

# Workaround for TeamCity CI, temp folders are created with root permissions, unacessible for USER
#RUN groupadd -r washu && useradd -ms /bin/bash -g washu user
#WORKDIR /home/user
# USER user

# Download Picard tools
RUN cd ~ && wget -q https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar

# Download ZINBRA
RUN cd ~ && wget -q https://github.com/JetBrains-Research/zinbra/releases/download/v0.4.0/zinbra-0.4.0.jar -O zinbra.jar
# Alternative CI link
# https://teamcity.jetbrains.com/repository/download/Epigenome_Zinbra/lastPinned/zinbra-0.4.0.jar?guest=1

