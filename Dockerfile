# DSS Methylation Analysis Container
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV R_VERSION=4.3.2
ENV PYTHON_VERSION=3.11

# Update system and install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    software-properties-common \
    apt-transport-https \
    ca-certificates \
    gnupg \
    lsb-release \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libhdf5-dev \
    libgsl-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

# Install R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y r-base r-base-dev r-base-core r-recommended && \
    rm -rf /var/lib/apt/lists/*

# Install additional R development libraries needed for rpy2
RUN apt-get update && apt-get install -y \
    libreadline-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python and pip
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Create symlinks for python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install setuptools and wheel first
RUN pip3 install --no-cache-dir --upgrade setuptools wheel

# Install JupyterLab and related packages
RUN pip3 install --no-cache-dir \
    jupyterlab==4.0.9 \
    jupyter-server-proxy \
    ipykernel \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    plotly

# Install conda 
RUN wget -qO- https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/miniforge
ENV PATH="/opt/miniforge/bin:$PATH"

# Initialize conda
RUN conda init bash

# Install bioconda packages for long-read methylation analysis
RUN conda install -y -c conda-forge -c bioconda \
    samtools \
    bcftools \
    htslib \
    bedtools \
    pysam

# Install tools via conda (use latest available versions)
RUN conda install -y -c bioconda pb-cpg-tools

# Install additional long-read analysis tools
RUN conda install -y -c bioconda \
    minimap2 \
    nanopolish



# Create symlinks for commonly used tools
RUN ln -s /opt/miniforge/bin/pb-cpg-tools /usr/local/bin/pb-cpg-tools && \
    ln -s /opt/miniforge/bin/samtools /usr/local/bin/samtools && \
    ln -s /opt/miniforge/bin/bedtools /usr/local/bin/bedtools


# Note: rpy2 removed due to compatibility issues - R and Python can be used separately

# Install R packages for methylation analysis
RUN R -e "install.packages(c('BiocManager', 'devtools', 'IRkernel'), repos='https://cloud.r-project.org')"

# Install Bioconductor and DSS
RUN R -e "BiocManager::install(c('DSS', 'bsseq', 'GenomicRanges', 'IRanges', 'methylKit'))"

# Install additional useful R packages
RUN R -e "install.packages(c('data.table', 'ggplot2', 'dplyr', 'readr', 'tidyverse', 'reshape2', 'corrplot', 'pheatmap', 'RColorBrewer', 'gridExtra'), repos='https://cloud.r-project.org')"

# Install R kernel for Jupyter
RUN R -e "IRkernel::installspec(user = FALSE)"

RUN conda install -y -c bioconda ont-modkit

# Create working directory
WORKDIR /workspace

# Create directories for data and results
RUN mkdir -p /workspace/data /workspace/results /workspace/scripts

# Set up JupyterLab configuration
RUN jupyter lab --generate-config

# Create a startup script
RUN echo '#!/bin/bash\n\
echo "Starting DSS Methylation Analysis Environment"\n\
echo "======================================"\n\
echo "Available kernels:"\n\
jupyter kernelspec list\n\
echo "======================================"\n\
echo "Long-read methylation tools:"\n\
echo "- modkit: $(modkit --version 2>/dev/null || echo "installed")"\n\
echo "- pb-cpg-tools: $(pb-cpg-tools --version 2>/dev/null || echo "installed")"\n\
echo "- samtools: $(samtools --version | head -1)"\n\
echo "- bedtools: $(bedtools --version)"\n\
echo "- minimap2: $(minimap2 --version 2>/dev/null || echo "installed")"\n\
echo "======================================"\n\
echo "Access JupyterLab at: http://localhost:8888"\n\
echo "======================================"\n\
exec jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --NotebookApp.token="" --NotebookApp.password=""' > /start.sh

RUN chmod +x /start.sh

# Expose JupyterLab port
EXPOSE 8888

# Set the default command
CMD ["/start.sh"]