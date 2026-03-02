# syntax=docker/dockerfile:1
# ============================================
# Stage 1: Build Angular frontend
# ============================================
FROM node:22-alpine AS frontend-build

WORKDIR /app/frontend

# Copy package files and install dependencies
COPY web/package*.json ./
RUN --mount=type=cache,target=/root/.npm \
    npm ci

# Copy frontend source and build
COPY web/ ./
RUN npm run build -- --configuration production

# ============================================
# Stage 2: Builder - Install Python dependencies
# ============================================
FROM python:3.10-slim AS backend-builder

# Prevent prompts
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1

# Install build dependencies (only needed for compiling)
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libffi-dev \
    wget

# Create virtual environment to isolate dependencies
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade core build tools
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip setuptools wheel "Cython<3.0"

# Copy and install Python dependencies
COPY backend/requirements.txt .
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install -r requirements.txt uvicorn

# Download MUSCLE5 (linux_intel64 runs on both amd64 and arm64 via emulation)
RUN wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64 -O /tmp/muscle \
    && chmod +x /tmp/muscle

# ============================================
# Stage 3: Runtime - Final image with everything
# ============================================
FROM python:3.10-slim

# Runtime environment variables
ENV PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH"

# Install runtime dependencies including nginx and supervisor
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libgomp1 \
    nginx \
    supervisor

# Copy virtual environment from backend-builder
COPY --from=backend-builder /opt/venv /opt/venv

# Copy MUSCLE binary from backend-builder
COPY --from=backend-builder /tmp/muscle /usr/local/bin/muscle

# Set working directory
WORKDIR /app

# Copy backend application code
COPY backend/ ./backend/

# Copy built Angular app from frontend-build stage
COPY --from=frontend-build /app/frontend/dist/web/browser /usr/share/nginx/html

# Copy nginx configuration
COPY nginx.conf /etc/nginx/nginx.conf

# Copy supervisor configuration
COPY supervisord.conf /etc/supervisor/conf.d/supervisord.conf

# Create necessary directories
RUN mkdir -p /var/log/supervisor

# Expose port 80 (nginx will handle both frontend and API proxy)
EXPOSE 80

# Use supervisor to manage both nginx and FastAPI
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]
