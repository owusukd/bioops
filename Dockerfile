# Multi-stage build for BioOps
FROM ubuntu:22.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install system dependencies and build tools
RUN apt-get update && apt-get install -y \
    build-essential \
    dpkg-dev \
    unzip \
    cmake \
    git \
    tar \
    libboost-dev \
    libarpack2-dev \
    libf2c2-dev \
    libeigen3-dev \
    python3-all \
    python3-all-dev \
    python3-pip \
    libopenblas-serial-dev \
    liblapack-dev \
    libsuitesparse-dev \
    libsuperlu-dev \
    python3-venv \
    libssl-dev \
    libffi-dev \
    curl \
    openjdk-17-jdk \
    openjdk-17-jre \
    wget \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Set Java environment
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# Create tools directory
WORKDIR /tools

# Install FATCAT
RUN git clone https://github.com/GodzikLab/FATCAT-dist.git && \
    chmod +x FATCAT-dist/FATCATMain/*

# Install P2Rank
RUN wget https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz && \
    tar -xzf p2rank_2.5.tar.gz && \
    rm p2rank_2.5.tar.gz && \
    chmod +x p2rank_2.5/prank

# Install APBS
RUN wget https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Linux.zip && \
    unzip APBS-3.4.1.Linux.zip && \
    rm APBS-3.4.1.Linux.zip && \
    chmod +x APBS-3.4.1.Linux/bin/*

# Install APoc (assuming it's available in the repo)
COPY apoc_v1b16.tar.gz /tools/
RUN tar -xzf apoc_v1b16.tar.gz && \
    rm apoc_v1b16.tar.gz && \
    chmod +x apoc/bin/*

# Install PDB2PQR
RUN apt-get update && apt-get install -y pdb2pqr && \
    rm -rf /var/lib/apt/lists/*

# Final stage
FROM ubuntu:22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    openjdk-17-jre \
    libboost-all-dev \
    libarpack2 \
    libf2c2 \
    libopenblas0 \
    liblapack3 \
    libsuperlu5 \
    pdb2pqr \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set Java environment
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# Copy tools from builder
COPY --from=builder /tools /tools

# Set up tool paths
ENV PATH="/tools/FATCAT-dist/FATCATMain:/tools/p2rank_2.5:/tools/APBS-3.4.1.Linux/bin:/tools/apoc/bin:$PATH"

# Create app directory
WORKDIR /app

# Create necessary directories
RUN mkdir -p logs tempDownloadDir/results config

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy logging files
COPY logging_config.py ./

# Copy application files
COPY bioops.py protFuncs.py ./
COPY logo/ ./logo/
COPY config/ ./config/

# Create non-root user
RUN useradd -m -u 1000 bioops && \
    chown -R bioops:bioops /app

# Switch to non-root user
USER bioops

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Set Streamlit configuration
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

# Run the application
CMD ["streamlit", "run", "bioops.py", "--server.port=8501", "--server.address=0.0.0.0"]
