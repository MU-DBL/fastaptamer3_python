# Use official lightweight Python image
FROM --platform=linux/amd64 python:3.10-slim

# Prevent prompts; speed up pip
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1

# OS build/runtime libs needed
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    zlib1g-dev libbz2-dev liblzma-dev libffi-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install MUSCLE5
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64 -O /usr/local/bin/muscle \
    && chmod +x /usr/local/bin/muscle

# Upgrade core build tools
RUN pip install --upgrade pip setuptools wheel "Cython<3.0"

# Set working directory
WORKDIR /app

# Copy dependency file first for better caching
COPY requirements.txt .

# Install dependencies
RUN pip install -r requirements.txt

# Copy app source code
COPY . .

# FastAPI port
EXPOSE 5001

# Start FastAPI
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "5001"]