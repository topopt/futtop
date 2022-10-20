SHELL   = /bin/bash
CC     ?= gcc
FUTC   ?= futhark
CFLAGS ?= -std=c99 -fopenmp -O3 -march=native

C_LIBS    = -std=c99 -lm
MC_LIBS   = -std=c99 -lm -lpthread
OCL_LIBS  = -lm -lOpenCL
CUDA_LIBS = -lm -lcuda -lcudart -lnvrtc

all: futtop_c

futtop: futtop.c libmultigrid.o io.o
	$(CC) $(CFLAGS) -o $@ $^ $(MC_LIBS)

libmultigrid.c: libmultigrid.fut src/*.fut
	$(FUTC) multicore --library $<

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $$(ls src/*.fut | cut -f1 -d'.')
	rm -f futtop futtop_library.c futtop_library.h
	rm -f src/*.c

bench:
	$(FUTC) bench --backend=opencl src/*.fut
