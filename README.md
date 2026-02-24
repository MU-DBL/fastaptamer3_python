# Build the image
docker build --no-cache -t yongfangqin/fastaptamer3:1.0 .
docker run -d -p 80:80 --name fastaptamer3 yongfangqin/fastaptamer3:1.0


# Build and push multi-platform image (run once to publish)
docker buildx create --use
docker buildx build -f Dockerfile --platform linux/amd64,linux/arm64 -t yongfangqin/fastaptamer3:1.1 --push .


# Package to tar
docker save -o fastaptamer3-frontend.tar yongfangqin/fastaptamer3-frontend:latest
docker save -o fastaptamer3-backend.tar yongfangqin/fastaptamer3-backend:latest
