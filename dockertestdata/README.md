Docker Image with test data
==========

This is just an image `continuumio/miniconda` with test data from
`/mnt/stripe/washu_test_data`.

Build
-------
To build you have to copy file `washu_test_data.tar.gz` to folder with this image.

```bash
# Copy test data to folder with Docker file
tar -cvzf washu_test_data.tar.gz -C /mnt/stripe/ washu_test_data

# Build Docker
docker build -t biolabs/test-data .

# Clean up
rm washu_test_data.tar.gz
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