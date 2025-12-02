# Build image
docker build -f Dockerfile -t fastaptamer3 .

# Run container
docker run -d -p 5001:5001 fastaptamer3
