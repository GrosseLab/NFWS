FROM broadinstitute/genomes-in-the-cloud:2.2.3-1469027018
MAINTAINER Martin Porsch <martin.porsch@informatik.uni-halle.de>

# install system tools
RUN apt-get update -qq
RUN apt-get upgrade -qqy
RUN apt-get install -qqy git pigz
RUN apt-get clean

# install conda, python, snakemake, rpy2 and R
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh &&  bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
ENV PATH /opt/conda/bin:$PATH
RUN conda create -n snakemake python=3.4 sqlite && /bin/bash -c "source activate snakemake" && conda install -c johanneskoester -c r snakemake rpy2

# install R-Packages
RUN R -e 'install.packages(c("ggplot2", "data.table"), repos="http://cran.cnr.Berkeley.edu")'

# pulling further install scripts and running them
RUN apt-get install make
RUN git clone https://github.com/GrosseLab/InstallProcedures.git #1
RUN InstallProcedures/bwa_0.7.15.sh
RUN InstallProcedures/bedtools_2.22.1.sh
RUN InstallProcedures/bedtools_2.25.sh
RUN InstallProcedures/samtools_1.3.1.sh
RUN InstallProcedures/fastqc_0.11.3.sh
RUN InstallProcedures/subread_1.5.0-p3.sh
RUN InstallProcedures/picard_2.5.0.sh
RUN cp InstallProcedures/picard /usr/local/bin/
RUN mkdir /usr/local/bin/GenomeAnalysisTK-3.6 && mv /usr/gitc/GATK36.jar /usr/local/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar && cp InstallProcedures/gatk /usr/local/bin/

# pulling core-facility rules
RUN mkdir /data && cd /data && git clone --branch develop https://github.com/GrosseLab/DefaultPipelines.git

# cleanup
#RUN apt-get autoremove -qqy python3 r-base r-base-dev
RUN rm -rf /usr/gitc

# startup
VOLUME ["/data/in"]
WORKDIR /data/in
ENTRYPOINT ["snakemake", "-j", "-k", "-p"]
