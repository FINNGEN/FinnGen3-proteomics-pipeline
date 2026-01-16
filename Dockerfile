# Dockerfile for FinnGen 3 Olink Analysis Pipeline
# Refactored version - Docker-ready for public release
# Version: 1.2.0

FROM rocker/r-ver:4.3.2

# Add version label for tracking
LABEL version="1.2.0" \
      maintainer="Reza Jabal <rjabal@broadinstitute.org>" \
      description="FinnGen 3 Olink Proteomics Analysis Pipeline"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    libv8-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    python3 \
    python3-pip \
    git \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install Google Cloud SDK for GCS access (required for pQTL steps)
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    rm -rf /var/lib/apt/lists/*

# Install PLINK (required for pQTL steps)
# Use explicit error handling and verify each step
RUN set -eux && \
    wget --no-verbose --show-progress https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip && \
    unzip -q plink_linux_x86_64_20230116.zip && \
    mv plink /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink && \
    rm -f plink_linux_x86_64_20230116.zip plink.log && \
    plink --version

# Set working directory
WORKDIR /pipeline

# Copy package installation script
COPY install_packages.R /pipeline/

# Install R packages
RUN Rscript install_packages.R

# Install BiocManager and BioConductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org/'); \
          BiocManager::install('sva', update=FALSE, ask=FALSE)"

# Copy pipeline scripts
COPY scripts/ /pipeline/scripts/
COPY config/ /pipeline/config/

# Create output and data directories
RUN mkdir -p /pipeline/output /pipeline/data

# Set environment variables
ENV PIPELINE_HOME=/pipeline
ENV PATH=/usr/local/bin:$PATH

# Default command
CMD ["Rscript", "/pipeline/scripts/run_pipeline.R", "--help"]




