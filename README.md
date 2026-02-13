# Build the image
docker build --no-cache -t yongfangqin/fastaptamer3:latest .
docker run -d -p 80:80 --name fastaptamer3 yongfangqin/fastaptamer3:latest


# Package to tar
docker save -o fastaptamer3-frontend.tar yongfangqin/fastaptamer3-frontend:latest
docker save -o fastaptamer3-backend.tar yongfangqin/fastaptamer3-backend:latest
