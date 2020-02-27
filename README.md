* How to build

create a bin subdirectory: i.e. by typing `mkdir build`

install an libSDL2, libpng, mpi and OpenCL on your distribution together with
the development version of packages (i.e. with header files)

MPI has been tested to work with OpenMPI but may work with other version.

OpenCL has been tested with the version shipped with the Nvidia Cuda Toolkit 
and the Intel OpenCL SDK. The later is only available for code generation on 
CPUs and hence only interesting for debugging purposes where an Nvidia GPU is 
not available.

edit the makefile, i.e. `emacs makefile`
adjust the header of the makefile to point to the correct locations 
of the libraries on your distribitution.

type `make`

if everything goes well you should have all the tools available in the bin
folder.