# SAIGE-QTL GPU Docker Container - GitHub Version
# This version clones code from GitHub during build
# Use this when your code is pushed to GitHub (recommended)

# Stage 1: Build environment
FROM nvidia/cuda:12.6.0-devel-ubuntu22.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libpcre2-dev \
    libreadline-dev \
    tzdata \
    libtbb-dev \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R 4.4.0
RUN wget https://cran.r-project.org/src/base/R-4/R-4.4.0.tar.gz && \
    tar -xzf R-4.4.0.tar.gz && \
    cd R-4.4.0 && \
    ./configure --prefix=/usr/local \
                --enable-R-shlib \
                --with-blas \
                --with-lapack \
                --with-readline=yes \
                --with-x=no && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf R-4.4.0 R-4.4.0.tar.gz

# Set CUDA environment variables
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=$CUDA_HOME/bin:$PATH
ENV LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
ENV CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME

# Clone repo first to get install_packages.R
WORKDIR /opt
RUN git clone https://github.com/exascale-genomics/SAIGEQTL-GPU.git

# Set MPI configuration for OpenMPI (not MPICH like on Polaris)
ENV MPI_TYPE=OPENMPI
ENV MPICH_GPU_SUPPORT_ENABLED=0

# Run the install script
WORKDIR /opt/SAIGEQTL-GPU
RUN Rscript ./extdata/install_packages.R

# Install pbdMPI with OpenMPI
RUN R -e "install.packages('pbdMPI', \
    configure.args='--with-mpi-type=OPENMPI', \
    repos='https://cloud.r-project.org/')"

# Build the package
RUN R CMD INSTALL --build .

# Stage 2: Runtime environment (smaller final image)
FROM nvidia/cuda:12.6.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# Install runtime dependencies only
RUN apt-get update && apt-get install -y \
    libopenblas0 \
    libgomp1 \
    libgfortran5 \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    libssl3 \
    libxml2 \
    zlib1g \
    libreadline8 \
    tzdata \
    libtbb-dev \
    openmpi-bin \
    libopenmpi3 \
    && rm -rf /var/lib/apt/lists/*

# Copy R installation from builder
COPY --from=builder /usr/local /usr/local

# Copy CUDA libraries needed for runtime
COPY --from=builder /usr/local/cuda/lib64 /usr/local/cuda/lib64

# Copy SAIGE-QTL installation
COPY --from=builder /opt/SAIGE-QTL /opt/SAIGE-QTL
COPY --from=builder /usr/local/lib/R /usr/local/lib/R

# Set environment variables
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=/opt/SAIGE-QTL/SAIGEQTL/extdata:$CUDA_HOME/bin:$PATH
ENV LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
ENV R_LIBS=/usr/local/lib/R/site-library:/usr/local/lib/R/library
ENV MPICH_GPU_SUPPORT_ENABLED=0

# Create working directory
WORKDIR /data

# Set entrypoint
ENTRYPOINT ["/bin/bash"]

# Labels
LABEL maintainer="Alex Rodriguez"
LABEL description="GPU-accelerated SAIGE-QTL for eQTL analysis with repeated measures"
LABEL version="1.0"
LABEL source="https://github.com/exascale-genomics/SAIGEQTL-GPU.git"
