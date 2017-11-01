Docker Image with test data
==========

This is just an image `continuumio/miniconda` with test data from
`/mnt/stripe/chip-seq-pipeline-test-data`.

Build
-------
To build you have to copy file to folder with image.

```bash
# Copy test data to folder with Docker file
tar -cvzf chip-seq-pipeline-test-data.tar.gz -C /mnt/stripe/ chip-seq-pipeline-test-data

# Build Docker
docker build -t biolabs/test-data .

# Clean up
rm chip-seq-pipeline-test-data.tar.gz
```

Push
-------
Before push you have to login to docker hub first.
```bash
docker login -u biolabs
```

Then you just push current image 
```bash
docker push biolabs/test-data
```