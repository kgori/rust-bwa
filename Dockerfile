FROM rust:bookworm

WORKDIR /usr/src/rust-bwa

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libclang-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . .

RUN cargo build --release
