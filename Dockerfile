# SAIGE-QTL GPU Docker Container
# Multi-stage build for efficient image size

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
    pkg-config \
    libopenblas-dev \
    liblapack-dev \
    libbz2-dev \
    libzstd-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libpcre2-dev \
    libreadline-dev \
    libpng-dev \
    libgit2-dev \
    libpq-dev \
    libmariadb-dev \
    tzdata \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    libsuperlu-dev \
    && rm -rf /var/lib/apt/lists/*

# Install older TBB version (2020.3) that has concurrent_vector
WORKDIR /tmp
RUN wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2020.3.tar.gz && \
    tar -xzf v2020.3.tar.gz && \
    cd oneTBB-2020.3 && \
    make -j$(nproc) compiler=gcc && \
    mkdir -p /usr/local/include/tbb && \
    cp -r include/tbb/* /usr/local/include/tbb/ && \
    find build -name "*.so*" -exec cp {} /usr/local/lib/ \; && \
    ldconfig && \
    cd /tmp && \
    rm -rf oneTBB-2020.3 v2020.3.tar.gz

# Install R 4.4.0
WORKDIR /tmp
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
    cd /tmp && \
    rm -rf R-4.4.0 R-4.4.0.tar.gz

# Set CUDA environment variables
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=$CUDA_HOME/bin:$PATH
ENV LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
ENV CUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME

# Set MPI configuration for OpenMPI
ENV MPI_TYPE=OPENMPI
ENV MPICH_GPU_SUPPORT_ENABLED=0

# Clone SAIGE-QTL from GitHub
WORKDIR /opt
RUN git clone https://github.com/exascale-genomics/SAIGEQTL-GPU.git

# Navigate into the repository
WORKDIR /opt/SAIGEQTL-GPU

# Run the install_packages.R script
RUN Rscript ./extdata/install_packages.R

# Install savvy library (VCF/BCF reader) and dependencies
RUN mkdir -p thirdParty/cget/include && \
    cd thirdParty && \
    # First install shrinkwrap (dependency of savvy)
    git clone https://github.com/jonathonl/shrinkwrap.git && \
    cd shrinkwrap && \
    mkdir -p build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    # Now install savvy
    git clone https://github.com/statgen/savvy.git && \
    cd savvy && \
    mkdir -p build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    cp -r ../include/savvy ../../cget/include/ && \
    cd ../../..

# Install pbdMPI with OpenMPI configuration
RUN R -e "install.packages('pbdMPI', \
    configure.args='--with-mpi-type=OPENMPI', \
    repos='https://cloud.r-project.org/')"

# Debug: Show original Makevars
RUN echo "=== ORIGINAL MAKEVARS ===" && cat src/Makevars

# Fix Makevars comprehensively
RUN sed -i 's|^LOCAL_HEADERS = .*|LOCAL_HEADERS =|g' src/Makevars && \
    sed -i 's|^LOCAL_LIBS = .*|LOCAL_LIBS =|g' src/Makevars && \
    sed -i 's|^CUDA_HOME = .*|CUDA_HOME = /usr/local/cuda|g' src/Makevars && \
    sed -i 's|-I $(LOCAL_HEADERS)||g' src/Makevars && \
    sed -i 's|-L$(LOCAL_LIBS)||g' src/Makevars && \
    sed -i 's|\$(TBBROOT)|/usr/local|g' src/Makevars && \
    sed -i 's|\.\./thirdParty|./thirdParty|g' src/Makevars && \
    sed -i 's|/usr/lib/aarch64-linux-gnu|/usr/lib/x86_64-linux-gnu|g' src/Makevars && \
    sed -i 's|MPI_CPPFLAGS = .*|MPI_CPPFLAGS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi|g' src/Makevars && \
    sed -i 's|MPI_LDFLAGS = .*|MPI_LDFLAGS = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi|g' src/Makevars

# Debug: Show updated Makevars
RUN echo "=== UPDATED MAKEVARS ===" && cat src/Makevars

# Fix source code to include TBB header properly
RUN if ! grep -q "#include.*concurrent_vector" src/GENO_null.hpp; then \
        sed -i '1i #include <tbb/concurrent_vector.h>' src/GENO_null.hpp && \
        echo "Added TBB concurrent_vector header"; \
    fi

# Set global R build environment for headers/libs (fixes MPI, CUDA, TBB include)
ENV PKG_CPPFLAGS="-I/usr/local/include -I/usr/local/include/tbb -I/usr/local/cuda/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi"
ENV PKG_CXXFLAGS="-I/usr/local/include -I/usr/local/include/tbb -I/usr/local/cuda/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi"
ENV PKG_LIBS="-L/usr/local/lib -ltbb"

# Build SAIGEQTL package
RUN R CMD INSTALL --build .

# Stage 2: Runtime environment
FROM nvidia/cuda:12.6.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# Install runtime dependencies
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

# Copy from builder
COPY --from=builder /usr/local /usr/local
COPY --from=builder /usr/local/cuda/lib64 /usr/local/cuda/lib64
COPY --from=builder /opt/SAIGEQTL-GPU /opt/SAIGEQTL-GPU
COPY --from=builder /usr/local/lib/R /usr/local/lib/R
COPY --from=builder /usr/local/lib/libtbb* /usr/local/lib/

# Update library cache
RUN ldconfig

# Set environment variables
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=/opt/SAIGEQTL-GPU/extdata:$CUDA_HOME/bin:$PATH
ENV LD_LIBRARY_PATH=$CUDA_HOME/lib64:/usr/local/lib:$LD_LIBRARY_PATH
ENV R_LIBS=/usr/local/lib/R/site-library:/usr/local/lib/R/library
ENV MPICH_GPU_SUPPORT_ENABLED=0

WORKDIR /data
ENTRYPOINT ["/bin/bash"]

# Labels
LABEL maintainer="Alex Rodriguez"
LABEL description="GPU-accelerated SAIGE-QTL for eQTL analysis"
LABEL version="1.0"
LABEL source="https://github.com/exascale-genomics/SAIGEQTL-GPU"

