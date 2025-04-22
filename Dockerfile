FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system packages and dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    build-essential \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    samtools \
    python3-pip \
    git \
    r-base \
    hisat2 \
    fastqc \
    htseq \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c(\
  'DESeq2', 'pheatmap', 'ggplot2', 'clusterProfiler', \
  'org.At.tair.db', 'BiocManager'), repos='http://cran.us.r-project.org')"

# Set default command
CMD ["/bin/bash"]


