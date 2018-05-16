Docker Image with test data and tools
=====================================

This is just an image `ubuntu:latest` - Ubuntu LTS with all the environment setup. 

Build
-----
Build image:
```bash
docker build -t biolabs/washu .
```

Push
----
Before push you have to login to docker hub first.
```bash
docker login -u biolabs
```

Then you just push current image 
```bash
docker push biolabs/washu
```