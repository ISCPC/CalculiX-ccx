CFLAGS = -Wall -O3 -DARCH="Linux" -fopenmp -DAURORA -DSBLAS -D_TIMELOG
FFLAGS = -Wall -O3
LDFLAGS = -llapack -lheterosolver_openmp -lblas_sequential -lsblas_openmp -fopenmp
#CFLAGS = -Wall -O3 -DARCH="Linux" -fopenmp -DAURORA -DSBLAS -D_TIMELOG -ftrace -DDISABLE_INLINE
#FFLAGS = -Wall -O3 -ftrace
#LDFLAGS = -llapack -lheterosolver_openmp -lblas_sequential -lsblas_openmp -fopenmp -ftrace

CC=ncc
FC=nfort

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

#CSRCS=arpackcs_sub.c aurora_hs.c libhs.c timelog.c
CSRCS=libccx.c cgsolver.c timelog.c
FSRCS=op.f

COBJS = $(CSRCS:.c=.o)
FOBJS = $(FSRCS:.f=.o)

all: libccx.so

libccx.so: $(COBJS) $(FOBJS)
	$(CC) -fpic -shared -o libccx.so $(COBJS) $(FOBJS) $(LDFLAGS)

clean:
	-rm libccx.so
	-rm *.o

