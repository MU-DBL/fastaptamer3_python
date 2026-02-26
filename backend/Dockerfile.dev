# ============================================
# Stage 1: Builder - Install dependencies
# ============================================
FROM python:3.10-slim as builder

# Prevent prompts; speed up pip
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1

# Install build dependencies (only needed for compiling)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libffi-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment to isolate dependencies
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade core build tools
RUN pip install --upgrade pip setuptools wheel "Cython<3.0"

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Download MUSCLE5
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64 -O /tmp/muscle \
    && chmod +x /tmp/muscle

# ============================================
# Stage 2: Runtime - Minimal final image
# ============================================
FROM python:3.10-slim

# Runtime environment variables
ENV PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH"

# Install ONLY runtime dependencies (no build tools)
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv

# Copy MUSCLE binary from builder
COPY --from=builder /tmp/muscle /usr/local/bin/muscle

# Set working directory
WORKDIR /app

# Copy application code
COPY . .

# FastAPI port
EXPOSE 5001

# Start FastAPI
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "5001"]