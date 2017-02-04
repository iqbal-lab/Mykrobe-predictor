## Install Mykrobe predictor with docker

Download the docker toolbox https://www.docker.com/products/docker-toolbox or docker engine https://docs.docker.com/engine/installation/. 

## Pull mykrobe predictor impage

	docker pull phelimb/mykrobe_predictor

## Run the image

	docker run phelimb/mykrobe_predictor mykrobe --help

## Run mykrobe predictor 

To run mykrobe predictor on data on your system you'll need to mount the directory where the data is located within your docker container. 

e.g. if your data is located within /Users/phelimb/Documents/data then you can run mykrobe like this:

	$ ls /Users/phelimb/Documents/data
	test.fastq.gz

	docker run -v /Users/phelimb/Documents/data:/data phelimb/mykrobe_predictor mykrobe predict sampleid tb -1 /data/tests.fastq.gz

You can pipe output locally

	docker run -v /Users/phelimb/Documents/data:/data phelimb/mykrobe_predictor mykrobe predict sampleid tb -1 /data/tests.fastq.gz > out.json
