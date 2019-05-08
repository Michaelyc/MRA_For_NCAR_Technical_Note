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

# Uncomment "-g" if the debug mode is desired
DEBUG=#-g

include=-I./include $(if $(mklIncPath),-I$(mklIncPath),) $(if $(dlibIncPath),-I$(dlibIncPath),) $(if $(armaIncPath),-I$(armaIncPath),)
LIBS=-lm -l$(mklLib) -l$(dlibLib) -l$(armaLib) $(if $(mklLibPath),-L$(mklLibPath),) $(if $(dlibLibPath),-L$(dlibLibPath),) $(if $(armaLibPath),-L$(armaLibPath),)

CXXFLAGS=$(include) $(DEBUG) -std=c++11 -O2

objDIR=./obj
srcDIR=./src

SRC=$(shell find $(srcDIR) -name '*.cpp')
OBJ=$(patsubst $(srcDIR)/%.cpp,$(objDIR)/%.o,$(SRC))

$(objDIR)/%.o: $(srcDIR)/%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

MRA: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(OBJ) MRA
