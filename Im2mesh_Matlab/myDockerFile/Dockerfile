# ==========================================
# Stage 1: The Builder (Heavy)
# ==========================================
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install all build dependencies
RUN apt-get update && apt-get -y install \
    libeigen3-dev \
    libgmp-dev \
    libgmpxx4ldbl \
    libmpfr-dev \
    git \
    make \
    cmake \
    gcc \
    g++

# Clone and build fTetWild
RUN git clone https://github.com/wildmeshing/fTetWild.git /fTetWild && \
    cd /fTetWild && \
    mkdir build && \
    cd build && \
    cmake .. && \
    # Using -j$(nproc) uses all CPU cores to speed up the build time!
    make -j$(nproc)

# ==========================================
# Stage 2: The Runner (Lightweight)
# ==========================================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Added libtbb12 here!
RUN apt-get update && apt-get -y install \
    libgmp10 \
    libgmpxx4ldbl \
    libmpfr6 \
    libgomp1 \
    libtbb12 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /fTetWild/build/FloatTetwild_bin /usr/local/bin/ftetwild

ENTRYPOINT ["ftetwild"]