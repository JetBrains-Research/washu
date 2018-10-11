FROM ubuntu:latest

USER root

# Update all the packages
RUN apt-get update --fix-missing

# Install conda
RUN apt-get install -y curl bzip2 ca-certificates
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
# SICER is python2 library, force it!
RUN sed -i 's#python#/opt/conda/envs/py2.7/bin/python#g' /opt/conda/envs/py2.7/share/sicer-1.1-3/SICER.sh
RUN ln -sf /opt/conda/envs/py2.7/share/sicer-1.1-3/SICER.sh /usr/local/bin/SICER.sh

# Python 3.5 environment with all the bioinformatics and python packages
RUN conda env create -q --name py3.5 --file /root/conda_envs/py3.5.environment.yml
RUN source activate py3.5 && pip install teamcity-messages pytest-codestyle==1.3.0 && source deactivate

# Bio utils available in path
RUN conda env create -q --name bio --file /root/conda_envs/bio.environment.yml
ENV PATH /opt/conda/envs/bio/bin/:$PATH

# Essential development tools and headers
RUN apt-get install --yes --fix-missing gawk build-essential libgl1-mesa-dev bc unzip quota libssl-dev \
    libboost-dev libssh2-1-dev libxml2-dev gfortran tar wget rsync
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

# Download Picard tools
RUN cd ~ && curl --location https://github.com/broadinstitute/picard/releases/download/2.10.7/picard.jar \
    --output picard.jar && chmod 755 ./picard.jar

# Install phantom peak qual tools
RUN cd ~ && curl --location https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/\
phantompeakqualtools/ccQualityControl.v.1.1.tar.gz --output phantompeakqualtools.tar.gz \
    && tar xvf phantompeakqualtools.tar.gz

# Download SPAN Peak Analyzer, see details: http://artyomovlab.wustl.edu/aging/span.html
RUN cd ~ && curl --location https://download.jetbrains.com/biolabs/span/span-0.7.0.4254.jar --output span.jar \
    && chmod 755 ./span.jar

# Install gosu; needed for correct user access permissions settings; see https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
RUN gpg --keyserver ha.pool.sks-keyservers.net --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4
RUN curl -o /usr/local/bin/gosu -SL "https://github.com/tianon/gosu/releases/download/1.4/gosu-$(dpkg --print-architecture)" \
    && curl -o /usr/local/bin/gosu.asc -SL "https://github.com/tianon/gosu/releases/download/1.4/gosu-$(dpkg --print-architecture).asc" \
    && gpg --verify /usr/local/bin/gosu.asc \
    && rm /usr/local/bin/gosu.asc \
    && chmod +x /usr/local/bin/gosu

COPY entrypoint.sh /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

