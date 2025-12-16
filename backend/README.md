# Build image
docker build -f Dockerfile -t yongfangqin/fastaptamer3-backend .

# Run container
docker run -d -p 5001:5001 fastaptamer3
