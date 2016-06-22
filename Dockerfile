# a docker image based on Ubuntu with snakemake installed
FROM broadinstitute/genomes-in-the-cloud:2.2.2-1466113830
MAINTAINER Martin Porsch <martin.porsch@informatik.uni-halle.de>
# install system tools
RUN apt-get update -qq
RUN apt-get upgrade -qqy
RUN apt-get install -qqy git pigz
# install R
RUN R -e 'install.packages(c("ggplot2", "data.table", "parallel"), repos="http://cran.cnr.Berkeley.edu")'
# install pyhthon
# RUN apt-get install -qqy python3 python3-pip
RUN pip3 install https://bitbucket.org/rpy2/rpy2/get/version_2.8.x.tar.gz && rm -rf /root/.cache
# install snakemake
RUN apt-get install -qqy python3-docutils python3-flask
RUN easy_install3 snakemake
# pulling further install scripts and running them
RUN git clone https://github.com/GrosseLab/InstallProcedures.git
RUN InstallProcedures/bwa_0.7.15.sh
RUN InstallProcedures/bedtools_2.22.1.sh
RUN InstallProcedures/bedtools_2.25.sh
RUN InstallProcedures/samtools_1.2.sh
RUN InstallProcedures/fastqc_0.11.3.sh
RUN InstallProcedures/subread_1.4.6.sh
RUN InstallProcedures/vcftools_0.1.14.sh
RUN mkdir /usr/local/bin/GenomeAnalysisTK-3.6 && mv /usr/gitc/GATK36.jar /usr/local/bin/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar && cp InstallProcedures/gatk /usr/local/bin/
RUN mkdir /usr/local/bin/Picard-1.1080 && mv /usr/gitc/picard.jar /usr/local/bin/Picard-1.1080/ && cp InstallProcedures/picard /usr/local/bin/

# pulling core-facility rules
RUN mkdir data && cd data && git clone --branch Docker https://github.com/GrosseLab/DefaultPipelines.git
# startup
ENTRYPOINT ["/bin/bash"]
