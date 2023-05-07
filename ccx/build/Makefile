CFLAGS = -w -O2 -std=gnu99 -fopenmp -DARCH="Linux" -DMKL_ILP64 -DPARDISO -DLONGLONG -DARPACK -DMATRIXSTORAGE -DUSE_MT=1
FFLAGS = -w -O2 -fopenmp -fallow-argument-mismatch -i8

CC=icc
FC=ifort

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.20.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

LIBS = \
       ../../../ARPACK/libarpack_INTEL.a \
       -lpthread -lm

ccx_2.20_MT: $(OCCXMAIN) ccx_2.20_MT.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.20.c; $(FC) -qopenmp -nofor-main -o $@ $(OCCXMAIN) ccx_2.20_MT.a $(LIBS) -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread

ccx_2.20_MT.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
