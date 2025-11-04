# Use official lightweight Python image
FROM python:3.10-slim

# Prevent prompts; speed up pip
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    # Prefer binary wheels; if none exist, fail instead of trying to build
    PIP_ONLY_BINARY=:all:

# OS build/runtime libs needed if any sdist slips through
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential python3-dev \
    zlib1g-dev libbz2-dev liblzma-dev libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# Upgrade core build tools so wheel resolution works
RUN pip install --upgrade pip setuptools wheel Cython

# Set working directory
WORKDIR /app

# Copy dependency file first for better caching
COPY requirements.txt .

# Install dependencies (cutadapt will use manylinux wheel)
RUN pip install -r requirements.txt

# Copy app source code
COPY . .

# FastAPI port
EXPOSE 5001

# Start FastAPI
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "5001"]