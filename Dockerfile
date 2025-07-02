# PyMaSC Development Docker Container
# Purpose: Linux environment testing for CI preparation with full dependencies

# Use x86_64 platform explicitly to avoid ARM64 x86intrin.h issues
FROM --platform=linux/amd64 python:3.8-slim

# Install system dependencies including all pysam requirements
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    git \
    libbz2-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libssl-dev \
    libncurses5-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /pymasc

# Copy only necessary files for dependency installation
COPY setup.py .
COPY PyMaSC/__init__.py PyMaSC/

# Install Python dependencies in optimal order
# First install build dependencies
RUN pip install --upgrade pip && \
    pip install numpy cython

# Install scientific libraries
RUN pip install scipy matplotlib

# Install PyMaSC dependencies 
RUN pip install "pysam>=0.22.0" "bx-python>=0.10.0"

# Install testing and development tools
RUN pip install pytest pytest-cov

# Create BitArray directory structure and copy source files only
RUN mkdir -p external/BitArray
COPY external/BitArray/*.c external/BitArray/*.h external/BitArray/Makefile external/BitArray/

# Build BitArray library in Linux environment with detailed verification
RUN cd external/BitArray && \
    rm -f libbitarr.a *.o && \
    make CC="gcc -fPIC" libbitarr.a && \
    ls -la libbitarr.a && \
    ar -tv libbitarr.a && \
    nm libbitarr.a | grep bit_array_set_bit && \
    echo "BitArray build successful with Linux symbols" && \
    cd ../..

# Copy the rest of the source code (excluding BitArray build artifacts)
COPY . .
# Remove any macOS-compiled BitArray artifacts to force clean Linux build
RUN rm -f external/BitArray/*.o external/BitArray/*.a

# Build and install PyMaSC
RUN python setup.py build_ext --inplace && \
    pip install -e .

# Verify installation with basic tests
RUN python -c "import PyMaSC; print('PyMaSC import successful')" && \
    pymasc --version

# Default command for comprehensive testing
CMD ["bash", "-c", "echo 'Running PyMaSC Docker Tests...' && \
     python -c 'import PyMaSC; print(\"PyMaSC import successful\")' && \
     pymasc --version && \
     python -m pytest tests/unit/ -v"]