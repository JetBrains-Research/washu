# Base Image based on Ubuntu x64
FROM biocontainers/biocontainers:latest

USER root
# Fix default shell for CONDA source activate command.
RUN ln -snf /bin/bash /bin/sh
RUN apt-get update
RUN apt-get install --yes bc

# java-jdk as dependecy
RUN conda install fastqc
RUN conda install bwa
RUN conda install bowtie
RUN conda install star
RUN conda install samtools
RUN conda install bedtools
RUN conda install deeptools
RUN conda install sra-tools
RUN conda install rseg
RUN conda install ucsc-bedgraphtobigwig
RUN conda install ucsc-bedclip
RUN pip install multiqc

RUN conda install --channel r r

# Install python, numpy and dependancies to build MACS
RUN apt-get install --yes build-essential git python python-numpy python-dev cython
# Clone MACS repository and checkout the requested tag
RUN cd /tmp && git clone https://github.com/taoliu/MACS.git && cd MACS && git checkout v2.0.9
# Compile and install MACS
RUN cd /tmp/MACS && python setup.py install

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
RUN source activate py3.5 && conda install pandas numpy && conda install --channel conda-forge matplotlib-venn

# Download Picard tools
RUN cd ~ && wget -q https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar

# Download ZINBRA
RUN cd ~ && wget -q https://github.com/JetBrains-Research/zinbra/releases/download/v0.4.0/zinbra-0.4.0.jar -O zinbra.jar

