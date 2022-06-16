This is the docker image for Comhet identification scripts.

## Instruction
1. Download docker image from [DockerHub](https://hub.docker.com/repository/docker/yirugong/ch_script)
```{bash}
docker push yirugong/ch_script:v1
```
2. Create a directory with input data (csv file)
3. Run docker image with customized data directory
```{bash}
docker run -it --rm -v <absolute data directory>:/data yirugong/ch_script:v1
# eg. ~/data2:/data
```
4. Obtain the output files in the same data directory
