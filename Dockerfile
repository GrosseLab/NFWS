# a docker image based on Ubuntu with snakemake installed
FROM broadinstitute/genomes-in-the-cloud:2.2.2-1466113830
MAINTAINER Martin Porsch <martin.porsch@informatik.uni-halle.de>

# install system tools
RUN apt-get update -qq
RUN apt-get upgrade -qqy
RUN apt-get install -qqy git pigz
# RUN apt-get autoremove python python3 r-base r-base-dev
RUN apt-get clean

# install conda, python, snakemake, rpy2 and R
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh &&  bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
ENV PATH /opt/conda/bin:$PATH
RUN conda create -n snakemake python=3.4 sqlite && /bin/bash -c "source activate snakemake" && conda install -c johanneskoester -c r snakemake rpy2

# install R-Packages
RUN R -e 'install.packages(c("ggplot2", "data.table"), repos="http://cran.cnr.Berkeley.edu")'

# pulling further install scripts and running them
RUN git clone https://github.com/GrosseLab/InstallProcedures.git #1
RUN InstallProcedures/bwa_0.7.15.sh
RUN InstallProcedures/bedtools_2.22.1.sh
RUN InstallProcedures/bedtools_2.25.sh
RUN InstallProcedures/samtools_1.2.sh
RUN InstallProcedures/fastqc_0.11.3.sh
RUN InstallProcedures/subread_1.5.0-p3.sh
RUN InstallProcedures/vcftools_0.1.14.sh
RUN mkdir /usr/local/bin/GenomeAnalysisTK-3.6 && mv /usr/gitc/GATK36.jar /usr/local/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar && cp InstallProcedures/gatk /usr/local/bin/
RUN mkdir /usr/local/bin/Picard-1.1080 && mv /usr/gitc/picard.jar /usr/local/bin/Picard-1.1080/ && cp InstallProcedures/picard /usr/local/bin/
RUN rm -rf /usr/gitc

# pulling core-facility rules
RUN mkdir /data && cd /data && git clone --branch develop https://github.com/GrosseLab/DefaultPipelines.git #

# startup
VOLUME ["/data/in"]
WORKDIR /data/in
ENTRYPOINT ["snakemake", "-j", "-k", "-p"]
