# Use Python 3.10 slim base image
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bedtools \
    git \
    build-essential \
    curl \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Create symlink for bedtools (epi2hit expects it at /usr/local/bin/bedtools)
RUN ln -sf /usr/bin/bedtools /usr/local/bin/bedtools

# Set working directory
WORKDIR /analysis

# Upgrade pip, setuptools, and wheel for better build support
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python packages
# Install core data science libraries first
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    scikit-learn

# Install visualization and analysis libraries
RUN pip install --no-cache-dir \
    matplotlib \
    seaborn \
    plotnine \
    statsmodels \
    statannot \
    joblib \
    tqdm \
    pyyaml

# Install packages that require compilation (with proper build support)
RUN pip install --no-cache-dir \
    pyBigWig \
    pybedtools

# Install epi2hit and biallelic_py
RUN pip install --no-cache-dir epi2hit==0.1.2

RUN pip install --no-cache-dir \
    git+https://bitbucket.org/weischenfeldt/biallelic_py.git@v0.1.0_epi2hit

# Copy entry point scripts
COPY epi_run_analysis.py /analysis/
COPY generate_metadata.py /analysis/
COPY entrypoint.sh /analysis/

# Make scripts executable
RUN chmod +x /analysis/entrypoint.sh /analysis/epi_run_analysis.py /analysis/generate_metadata.py

# Create data and results directories
RUN mkdir -p /analysis/data/results

# Set the default entry point
ENTRYPOINT ["/analysis/entrypoint.sh"]
CMD ["--help"]
