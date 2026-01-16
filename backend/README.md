# Build image
docker build -f Dockerfile --platform=linux/amd64 -t yongfangqin/fastaptamer3-backend .

# Run container
docker run -d -p 5001:5001 yongfangqin/fastaptamer3-backend
