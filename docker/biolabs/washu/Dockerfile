FROM ubuntu:latest

USER root

# Install conda
RUN apt-get update --fix-missing && apt-get install -y curl bzip2
RUN curl --location https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh --output ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH
RUN ln -snf /bin/bash /bin/sh

# Conda envs
COPY conda_envs /root/conda_envs

## Packages to be loaded using module in dedicated envs
RUN conda env create -q --name samtools --file /root/conda_envs/samtools.environment.yml
RUN conda env create -q --name bedtools --file /root/conda_envs/bedtools.environment.yml
RUN conda env create -q --name r --file /root/conda_envs/r.environment.yml
RUN conda env create -q --name bowtie --file /root/conda_envs/bowtie.environment.yml
RUN conda env create -q --name bowtie2 --file /root/conda_envs/bowtie2.environment.yml
RUN conda env create -q --name java --file /root/conda_envs/java.environment.yml
RUN conda env create -q --name sratoolkit --file /root/conda_envs/sratoolkit.environment.yml

# Python 2.7 with MACS2, SICER
RUN conda env create -q --name py2.7 --file /root/conda_envs/py2.7.environment.yml
RUN ln -sf /opt/conda/envs/py2.7/bin/macs2 /usr/local/bin/macs2
RUN ln -sf /opt/conda/envs/py2.7/bin/SICER.sh /usr/local/bin/SICER.sh

# Python 3.5 environment with all the bioinformatics and python packages
RUN conda env create -q --name py3.5 --file /root/conda_envs/py3.5.environment.yml
RUN source activate py3.5 && pip install teamcity-messages && source deactivate

# Bio utils available in path
RUN conda env create -q --name bio --file /root/conda_envs/bio.environment.yml
ENV PATH /opt/conda/envs/bio/bin/:$PATH

# Essential development tools and headers
RUN apt-get install --yes gawk build-essential libgl1-mesa-dev bc unzip quota libssl-dev \
    libboost-dev curl libssh2-1-dev libcurl4-openssl-dev libxml2-dev gfortran tar
# Fix missing gtar for R devtools
RUN ln -snf /bin/tar /bin/gtar

# Install spp R package for phantom peak qual tools
RUN source activate r && \
    R -e "source('https://bioconductor.org/biocLite.R');\
    biocLite(c('caTools', 'snow', 'snowfall', 'bitops', 'Rsamtools'));\
    install.packages('devtools', repos='http://cran.us.r-project.org');\
    require(devtools);options(unzip = 'internal');\
    devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)" && \
    source deactivate

# We need this for "which module" command
COPY ./module.sh /opt/
RUN ln -s /bin/echo /usr/bin/module

# Make new user and add it to the washu group
# TODO: uncomment when 'root' issue has been fixed in TC
#RUN groupadd -r washu && useradd -ms /bin/bash -g washu user
#WORKDIR /home/user
#
## Move testdata data
#RUN mv /root/washu_test_data.tar.gz /home/user && chmod a+r /home/user/washu_test_data.tar.gz
#
## Ordinary user
#USER user

# Download Picard tools
RUN cd ~ && curl --location https://github.com/broadinstitute/picard/releases/download/2.10.7/picard.jar \
    --output picard.jar

# Download ZINBRA
RUN cd ~ && curl --location https://github.com/JetBrains-Research/zinbra/releases/download/v0.5.1/zinbra-0.5.1.jar \
    --output zinbra.jar

# Install phantom peak qual tools
RUN cd ~ && curl --location https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/\
phantompeakqualtools/ccQualityControl.v.1.1.tar.gz --output phantompeakqualtools.tar.gz \
    && tar xvf phantompeakqualtools.tar.gz