# Manual

* Please note that there is full manual with detailed build and usage 
  instruction as well as a tutorial available at
  https://gitlab.in2p3.fr/mnhn-tools/mnhn-tree-tools-manual/-/blob/master/manual.pdf

# How to build (Debain 10 / Ubuntu 20.04  with standard non optimzed libraries)
```
sudo apt-get install git build-essential libpng-dev libsdl2-dev liblapack-dev libopenmpi-dev libpocl-dev ocl-icd-opencl-dev pocl-opencl-icd

git clone https://github.com/haschka/mnhn-tree-tools

cd mnhn-tree-tools
mkdir bin
make all
cd bin

# make MNHN-TREE-TOOLS available from any folder ( this is temporary, you
# may modify your .bashrc and similar files to make this permanent
export PATH=$PATH:$PWD
```

# How to build (generic and optimized)

* create a bin subdirectory: i.e. by typing `mkdir bin`

* install libSDL2, libpng, mpi and OpenCL on your distribution/unix-system
together with the development version of packages (i.e. with header files)

MPI has been tested to work with OpenMPI but may work with other version.

OpenCL has been tested with the version shipped with the Nvidia Cuda Toolkit 
and the Intel OpenCL SDK. The later is only available for code generation on 
CPUs and hence only interesting for debugging purposes where a GPU is 
not available.

* edit the makefile, i.e. `emacs makefile`
adjust the header of the makefile to point to the correct locations 
of the libraries on your distribitution.

* type `make all`

if everything goes well you should have all the tools available in the bin
folder.

* Flaws

Currently only works on 64 bit machines where size_t and pointers are 64bits
long.
