#############################################################################
# Makefile for building: Poisson
# Craig Lage - UC Davis
# 26-Feb-16
#############################################################################

BOOST_DIR     = /usr/local
HDF5_DIR     = /usr/local

####### Compiler, tools and options

CC            = gcc
CXX           = g++
CFLAGS        = -pipe -O2 -Wall -W
CXXFLAGS      = -pipe -O2 -Wall -W
INCPATH       =  -I.  -I$(HDF5_DIR)/include -I$(BOOST_DIR)/include
LINK          = g++
LFLAGS        =
LIBS          = -L$(BOOST_DIR)/lib -L$(HDF5_DIR)/lib -lm -lz -lhdf5 -lhdf5_cpp -lpthread

####### Output directory

DEL_FILE      = rm -f
OBJECTS_DIR   = ./

####### Files

SOURCES       = poisson.cpp \
		multigrid.cpp \
		array3d.cpp \
		array2d.cpp \
		polygon.cpp \
		fileio.cpp \
		hdf5write.cpp
OBJECTS       = poisson.o \
		multigrid.o \
		array3d.o \
		array2d.o \
		polygon.o \
		fileio.o \
		hdf5write.o
DESTDIR       =
TARGET        = Poisson

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)


clean:
	-$(DEL_FILE) $(OBJECTS) $(TARGET)

####### Compile

poisson.o: poisson.cpp poisson.h \
		multigrid.h \
		globals.h \
		fileio.h \
		hdf5write.h \
		array3d.h \
		array2d.h \
		polygon.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o poisson.o poisson.cpp

multigrid.o: multigrid.cpp multigrid.h \
		globals.h \
		fileio.h \
		hdf5write.h \
		array3d.h \
		array2d.h \
		polygon.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o multigrid.o multigrid.cpp

array3d.o: array3d.cpp array3d.h \
		globals.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o array3d.o array3d.cpp

array2d.o: array2d.cpp array2d.h \
		globals.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o array2d.o array2d.cpp

polygon.o: polygon.cpp polygon.h \
		globals.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o polygon.o polygon.cpp

fileio.o: fileio.cpp fileio.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fileio.o fileio.cpp

hdf5write.o: hdf5write.cpp hdf5write.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hdf5write.o hdf5write.cpp
