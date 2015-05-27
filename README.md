Mykrobe-predictor
=================

The Mykrobe-predictor can either be [built and compiled from source](#compilation) or [run using a ready-made Docker image](#docker-usage).

### Download ###

This repo uses submodules so to download it, type the following

git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git

Then run the following command to compile libraries. 

sh install.sh 

### Compilation ###

To compile for S.aureus, run the following command 

**make STAPH=1 predictor**

to create this executable 

bin/Mykrobe.predictor.staph

To compile for M.tuberculosis, run the following command 

**make TB=1 predictor**

to create this executable 

bin/Mykrobe.predictor.tb


[In order to install on Mac you may need to install gcc via homebrew or macports. Running  "gcc -v" 
should return something like:
"gcc version 4.9.2 (Homebrew gcc49 4.9.2_1)"
not 
"Apple LLVM version 6.0 (clang-600.0.57)"

### Usage ###

#### Standard usage ####

Mykrobe.predictor.staph --file some_file.fastq --install_dir /path/to/Mykrobe-predictor

Finally, there are GUI versions of Mykrobe-predictor for Windows and Mac OS X, which you can download from Releases

#### Docker usage ####

FIrstly, ensure the Docker daemon is installed. See https://docs.docker.com/installation/#installation for installation guides fo Linux, Mac OS and Windows.

The simplest method to use the Dockerised CLI Mykrobe-predictor is to use the [ready-made image hosted at Docker Hub](https://registry.hub.docker.com/u/jetstack/mykrobe-predictor/). 
Using ```docker run```, the latest image will be automatically downloaded and run in an separate, isolated container and will only
live for the duration of the computation.

```shell
docker run -v /location/to/Mykrobe-predictor/data:/mykrobe/data -it jetstack/mykrobe-predictor staph --file data/path/to/data/file
```

* The executable is selected by keyword - e.g. ```staph``` or ```tb```.

* The Docker ```-v``` switch is used to map a directory on the host (your host PC/laptop) to the relevant directory *inside* the container (i.e. 
host-directory : container-directory). You only need specify the host directory and this should be a full absolute path.

* The Mykrobe-predictor ```--file``` and ```--list``` switches take the same form but the path to the data directory should be prefixed with data/.
(e.g. data/staph/antibiotics/blaZ.fa).

### Paper, citation ###
We have a preprint of the paper describing Mykrobe predictor here:
http://biorxiv.org/content/early/2015/04/26/018564
Please cite us if you use Mykrobe predictor in a publication
