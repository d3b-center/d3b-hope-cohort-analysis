FROM --platform=linux/amd64 rocker/tidyverse:4.4.0
LABEL maintainer = "Jo Lynne Rokita (rokita@chop.edu)"
WORKDIR /rocker-build/

### Install apt-getable packages to start
#########################################
RUN apt-get update  
RUN apt-get install -y --no-install-recommends apt-utils dialog
RUN apt-get install -y --no-install-recommends \
  libxt6 \
  bzip2 \
  libmagick++-dev

# Install dev libraries and curl
RUN apt update && apt install -y zlib1g-dev \
	libncurses5-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	curl \
	cmake

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19')"

RUN R -e 'BiocManager::install(c( \
	"Biobase", \
	"BiocManager", \
    "broom", \
	"circlize", \
    "ComplexHeatmap", \
	"corrplot", \
	"cowplot", \
	"cutpointr", \
    "data.table", \
	"ggalluvial", \
	"ggpubr", \
	"ggthemes", \
	"ggstatsplot", \
	"ggforce", \
	"ggfortify", \
	"ggrepel", \
	"GenomicFeatures", \
    "GenomicRanges", \
	"GSVA", \
	"msigdbr", \
	"NOISeq", \
	"optparse", \
    "reshape2", \
	"R.utils", \
	"survival", \
	"survminer", \
	"sva", \
	"txdbmaker", \
	"WGCNA" \
    ))'

# add GitHub R packages
RUN R -e "remotes::install_github('d3b-center/annoFusedata', ref = '321bc4f6db6e9a21358f0d09297142f6029ac7aa', dependencies = TRUE)"
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '90d64f8fc50bee7060be577f180ae019a9bbbb84', dependencies = TRUE)"
RUN R -e "devtools::install_github('andymckenzie/DGCA', dependencies = TRUE)"

# Install pip3 and low-level python installation reqs
RUN apt-get update
RUN apt-get -y --no-install-recommends install \
    python3-pip  python3-dev
RUN ln -s /usr/bin/python3 /usr/bin/python  
RUN python3 -m pip install --upgrade pip

RUN pip3 install \
    "Cython==0.29.15" \
    "setuptools==46.3.0" \
    "six==1.14.0" \
    "wheel==0.34.2" \
    "numpy==1.24.3" \
    "pandas==2.0.1" \
    "matplotlib==3.7.1" \
    "scikit-learn==1.2.2" \
    "seaborn==0.12.2" \
    "rpy2==3.5.0" \
    "utils==1.0.1" 
    
# Required for mapping segments to genes
# Add bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && rm -f bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/* /usr/local/bin && \
    cd .. && rm -rf bedtools2

# Add bedops per the BEDOPS documentation
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.37/bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    tar -jxvf bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    rm -f bedops_linux_x86_64-v2.4.37.tar.bz2 && \
    mv bin/* /usr/local/bin
    
ADD Dockerfile .
    
