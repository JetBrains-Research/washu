Docker Image with test data and tools
=====================================

This is just an image `biolabs/test-data` with all the tools installed. 

Build
-----

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
