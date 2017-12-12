Docker Image with test data and tools
=====================================

This is just an image `biolabs/test-data` with all the tools installed. 

Build
-----
Ensure your `biolabs/test-data` image is up to date.
```bash
docker pull biolabs/test-data
```

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

Export conda environment
----
```bash
conda env export --name ENV_NAME
```

Run tests localy
---
```bash
docker run -v ~/work/washu:/washu -it biolabs/washu
source activate py3.5 && cd /washu && bash test_pipeline.sh
```