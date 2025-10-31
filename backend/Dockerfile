# Use official lightweight Python image
FROM python:3.10-slim

# Set working directory
WORKDIR /app

ENV UPLOAD_DIR=files
ENV PYTHONUNBUFFERED=1

# Copy dependency file
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy app source code
COPY . .

# Expose FastAPI port
EXPOSE 5001

# Start FastAPI with uvicorn
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "5001"]
