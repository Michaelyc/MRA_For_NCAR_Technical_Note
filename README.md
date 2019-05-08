## About

Authour: Huang Huang (Email: <hhuang0402@gmail.com>)

Date: May 5, 2019

Description: This package contains the C++ source codes for the serial and parallel implementations of the multi-resolution approximation in Huang et al. (2019). More information about the methodology of the multi-resolution approximation can be seen in Katzfuss (2017).

[1] Huang, H., L. R. Blake, and D. M. Hammerling (2019). Pushing the limit: a hybrid parallel implementation of the multi-resolution approximation for massive data. NCAR Technical Note (NCAR/TN-558+STR).

[2] Katzfuss, M. (2017). A multi-resolution approximation for massive spatial datasets. Journal of the American Statistical Association 112(517), 201–214.

Redistribution: subject to the conditions in the provided "LICENSE", redistribution is allowed when citing Huang et al. (2019).

## Code structure

Directory "serial\_MRA" contains all the codes for the serial implementation and directory "parallel\_MRA" contains all the codes for the parallel implementation. User can decide which binary program to build based on the available computing facility. In either "serial\_MRA" or "parallel\_MRA",

1. Directory "src" contains the source files.
2. Directory "include" contains the header files.
3. Directory "obj" is an empty folder to store the objects to be built.
4. File "user\_parameters" is used to set up the parameters for the multi-resolution approximation execution. The instructions for setting each parameter are provided directly in "user\_parameters" in the lines starting with a hash sign "#".
5. File "Makefile" is used to compiling the code. See more discussion in Section [Compiling](#Compiling) later.


## Prequisite libraries

The following libraries should be installed in advance, which will be required for the program.

1. Intel MKL. See https://software.intel.com/en-us/mkl for installation tutorial.
2. armadillo. See http://arma.sourceforge.net for installation tutorial.
3. dlib. See http://www.dlib.net for installation tutorial.

For the "parallel_MRA" implementation, the MPI with OpenMP library is also required. Options include Intel MPI, Open MPI, SGI MPT. User just choose one that is available.

All the library names, paths, and corresponding header file paths should be included in "Makefile" correctly.

## <a name="Compiling"></a> Compiling 

### Compiling for the serial program

Enter the directory "serial_MRA" and assign correct values to the variables in "Makefile" below.

```
# Choose the available C++ compilier
CXX=

# Set the header file path for the prequisite libraries
## path for mkl header files
mklIncPath=
## path for dlib header files
dlibIncPath= 
## path for armadillo header files
armaIncPath= 

# Set the library path for the prequisite libraries
## path for mkl library
mklLibPath=
## path for dlib library
dlibLibPath=
## path for armadillo library
armaLibPath=

# Change the library names if what are installed have different names
## mkl library name
mklLib=mkl
## dlib library name
dlibLib=dlib
## armadillo library name
armaLib=armadillo
```

Then, run 
`make`
in the terminal to compile the codes. After compiling finishes, a binary program "MRA" will be generated. 

### Compiling for the parallel program

Enter the directory "parallel_MRA" and assign correct values to the variables in "Makefile" below.

```
# Choose the available MPI C++ compilier
CXX=

# OpenMP setting
## Set the header file path for OpenMP
ompIncPath=
## Set the header file path for OpenMP
ompLibPath=
## Change the OpenMP library name if what is installed has a different name
ompLib=omp

# Set the header file path for the prequisite libraries
## path for mkl header files
mklIncPath=
## path for dlib header files
dlibIncPath= 
## path for armadillo header files
armaIncPath= 

# Set the library path for the prequisite libraries
## path for mkl library
mklLibPath=
## path for dlib library
dlibLibPath=
## path for armadillo library
armaLibPath=

# Change the library names if what are installed have different names
## mkl library name
mklLib=mkl
## dlib library name
dlibLib=dlib
## armadillo library name
armaLib=armadillo
```

Then, run 
`make`
in the terminal to compile the codes. After compiling finishes, a binary program "MRA" will be generated. 

## Run the program
Carefully read the comments in "user_parameters" and assign the desired value to each parameter. 

1. In directory "serial_MRA", run
`MRA`
to get results.

2. In directory "parallel_MRA", use the correct MPI library command to run the program. For instance,
`mpirun -n 2 MRA`,
where 2 is the number of MPI processes specified.

## Contact
If you have any questions or face any problems for this code, please do not hesitate to [contact me](mailto:hhuang0402@gmail.com). 