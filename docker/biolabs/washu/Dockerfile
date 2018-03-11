FROM biolabs/test-data

USER root
# Fix default shell for CONDA source activate command.
RUN ln -snf /bin/bash /bin/sh

# Essential development tools and headers
RUN apt-get install --yes build-essential libgl1-mesa-dev bc unzip quota libssl-dev \
libboost-dev curl libssh2-1-dev libcurl4-openssl-dev

# GNU AWK requered for proper scripts work
RUN apt-get install --yes gawk

# Conda envs:
COPY conda_envs /root/conda_envs

# Hack to enable MACS2 in another conda environment
RUN conda env create -q --name macs2 --file /root/conda_envs/macs2.environment.yml
RUN ln -sf /opt/conda/envs/macs2/bin/macs2 /usr/local/bin/macs2

# SICER
RUN conda env update -q --file /root/conda_envs/root.environment.yml
RUN cd /tmp && curl --location http://home.gwu.edu/~wpeng/SICER_V1.1.tgz --output SICER_V1.1.tgz && \
    tar xvf SICER_V1.1.tgz && mv SICER_V1.1 /opt
# Please refer to README for installation instructions, modify scripts, i.e.
RUN sed -i 's#/home/data/SICER1.1#/opt/SICER_V1.1#g' /opt/SICER_V1.1/SICER/SICER.sh
# SICER is python2 library, force it!
RUN sed -i 's#python#python2#g' /opt/SICER_V1.1/SICER/SICER.sh
RUN chmod a+x /opt/SICER_V1.1/SICER/SICER.sh
ENV PATH $PATH:/opt/SICER_V1.1/SICER

# Install packages to be loaded using module in dedicated envs
RUN conda env create -q --name  samtools --file /root/conda_envs/samtools.environment.yml
RUN conda env create -q --name  bedtools --file /root/conda_envs/bedtools.environment.yml
RUN conda env create -q --name  r --file /root/conda_envs/r.environment.yml
RUN conda env create -q --name  bowtie --file /root/conda_envs/bowtie.environment.yml
RUN conda env create -q --name  bowtie2 --file /root/conda_envs/bowtie2.environment.yml
RUN conda env create -q --name  java --file /root/conda_envs/java.environment.yml
RUN conda env create -q --name  sratoolkit --file /root/conda_envs/sratoolkit.environment.yml

# Deeptools installs samtools and lots of dependencies, isolate it
RUN conda env create -q --name  deeptools --file /root/conda_envs/deeptools.environment.yml
RUN ln -s /opt/conda/envs/deeptools/bin/deeptools /usr/local/bin/deeptools
RUN ln -s /opt/conda/envs/deeptools/bin/bamCoverage /usr/local/bin/bamCoverage

# Install env py3.5
RUN conda env create -q --name py3.5 --file /root/conda_envs/py3.5.environment.yml

# Install spp R package for phantom peak qual tools
RUN /opt/conda/envs/r/bin/R -e "source('https://bioconductor.org/biocLite.R');biocLite(c('caTools', 'Rsamtools'));\
install.packages('devtools', repos='http://cran.us.r-project.org');require(devtools);\
devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)"

# To prevent problems with Java interfierence just move execuble to emulate module
RUN mkdir /opt/fastqc && mv /opt/conda/envs/java/bin/fastqc /opt/fastqc

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