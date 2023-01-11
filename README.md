# Futtop
This repository contains an implementation of a topology  optimisation solver for linear  elastic compliance minimisation in three dimensions. The implementation is based on the Futhark language.

## Compilation
The presented code has been tested using the Futhark compiler version 0.21.11, and GCC 11.2 with nvptx. 

The provided Makefile compiles the code for the multicore backend in Futhark by default. 

To enable the OpenCL backend modify the following lines in the Makefile:
```bash
futtop: futtop.c libmultigrid.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(MC_LIBS)

libmultigrid.c: libmultigrid.fut src/*.fut
	$(FUTC) multicore --library $<
```

to
```bash
futtop: futtop.c libmultigrid.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(OCL_LIBS)

libmultigrid.c: libmultigrid.fut src/*.fut
	$(FUTC) opencl --library $<
```

## Running the Code
The default design problem is a 2x1x1 cantilever problem. To run 20 iterations of the code on a grid of 128 times 64 times 64 voxels use the following commands:
```bash
$ ./top3d -x 16 -y 8 -z 8
```

A list of available options is printed on start of the program.

## Authorship
This code has been developed by Erik A. Träff under the supervision of Niels Aage and Ole Signmund.
