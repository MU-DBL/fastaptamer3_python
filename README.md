# Performance setup (run once before first use)
# Gives Docker more CPU/RAM for large SELEX file processing (~5 min → ~1 min).
# No admin rights required.

# Windows:
powershell -ExecutionPolicy Bypass -File setup-docker-resources.ps1
# Then restart Docker Desktop.

# macOS:
bash setup-docker-resources.sh
# Then restart Docker Desktop (menu bar icon → Restart).

# Build and push multi-platform image (run once to publish)
docker buildx create --use
docker buildx build -f Dockerfile --platform linux/amd64,linux/arm64 -t yongfangqin/fastaptamer3:1.1 --push .
docker buildx build --platform linux/amd64 -t yongfangqin/fastaptamer3:1.1 --load .

# Package to tar
docker save -o fastaptamer3-frontend.tar yongfangqin/fastaptamer3-frontend:latest
docker save -o fastaptamer3-backend.tar yongfangqin/fastaptamer3-backend:latest
