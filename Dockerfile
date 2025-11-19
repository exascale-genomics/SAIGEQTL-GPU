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
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Install older TBB version (2020.3) that has concurrent_vector
WORKDIR /tmp
RUN wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2020.3.tar.gz && \
    tar -xzf v2020.3.tar.gz && \
    cd oneTBB-2020.3 && \
    make -j$(nproc) && \
    cp -r include/tbb /usr/local/include/ && \
    cp build/*_release/*.so* /usr/local/lib/ && \
    ldconfig && \
    cd .. && \
    rm -rf oneTBB-2020.3 v2020.3.tar.gz

# Verify TBB installation
RUN ls -la /usr/local/include/tbb/ && \
    ls -la /usr/local/lib/libtbb*

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

# Set MPI configuration for OpenMPI
ENV MPI_TYPE=OPENMPI
ENV MPICH_GPU_SUPPORT_ENABLED=0

# Clone repo first to get install_packages.R
WORKDIR /opt
RUN git clone https://github.com/exascale-genomics/SAIGEQTL-GPU.git

# Run the install script
WORKDIR /opt/SAIGEQTL-GPU
RUN Rscript ./extdata/install_packages.R

# Install pbdMPI with OpenMPI
RUN R -e "install.packages('pbdMPI', \
    configure.args='--with-mpi-type=OPENMPI', \
    repos='https://cloud.r-project.org/')"
    
# Debug: Check if Makevars exists and show its content
RUN echo "=== Checking for Makevars ===" && \
    find . -name Makevars -type f && \
    if [ -f src/Makevars ]; then \
        echo "=== Original Makevars content ===" && \
        cat src/Makevars; \
    else \
        echo "=== No Makevars found, checking Makevars.in ===" && \
        if [ -f src/Makevars.in ]; then cat src/Makevars.in; fi; \
    fi

# Fix Makevars or Makevars.in to use correct TBB path
RUN if [ -f src/Makevars ]; then \
        sed -i 's|-I/usr/include/tbb|-I/usr/local/include|g' src/Makevars && \
        sed -i 's|-L/usr/lib/x86_64-linux-gnu|-L/usr/local/lib|g' src/Makevars && \
        echo "=== Updated Makevars ===" && cat src/Makevars; \
    fi && \
    if [ -f src/Makevars.in ]; then \
        sed -i 's|-I/usr/include/tbb|-I/usr/local/include|g' src/Makevars.in && \
        sed -i 's|-L/usr/lib/x86_64-linux-gnu|-L/usr/local/lib|g' src/Makevars.in && \
        echo "=== Updated Makevars.in ===" && cat src/Makevars.in; \
    fi

# Set R build environment variables
ENV PKG_CPPFLAGS="-I/usr/local/include -I/usr/local/include/tbb"
ENV PKG_CXXFLAGS="-I/usr/local/include -I/usr/local/include/tbb"
ENV PKG_LIBS="-L/usr/local/lib -ltbb"

# Build SAIGEQTL package
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
    openmpi-bin \
    libopenmpi3 \
    && rm -rf /var/lib/apt/lists/*

# Copy R installation from builder
COPY --from=builder /usr/local /usr/local

# Copy CUDA libraries needed for runtime
COPY --from=builder /usr/local/cuda/lib64 /usr/local/cuda/lib64

# Copy SAIGE-QTL installation (correct path: SAIGEQTL-GPU)
COPY --from=builder /opt/SAIGEQTL-GPU /opt/SAIGEQTL-GPU
COPY --from=builder /usr/local/lib/R /usr/local/lib/R

# Copy TBB libraries from builder
COPY --from=builder /usr/local/lib/libtbb* /usr/local/lib/

# Run ldconfig to update library cache
RUN ldconfig

# Set environment variables (correct paths)
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=/opt/SAIGEQTL-GPU/extdata:$CUDA_HOME/bin:$PATH
ENV LD_LIBRARY_PATH=$CUDA_HOME/lib64:/usr/local/lib:$LD_LIBRARY_PATH
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
