CC=mpicc
LD=${CC}
COVFLAGS = -fprofile-arcs -ftest-coverage
INCLUDES = -I$(OPENBLAS_ROOT)/include
CFLAGS+=-Wall -pedantic -O3 -fPIC $(INCLUDES) ${COVFLAGS}
DBGFLAGS = -g -pg

LDFLAGS+=-lm -fopenmp -L$(OPENBLAS_ROOT)/lib ${COVFLAGS} -g -pg -lopenblas

#OBJS=blas.o cg.o util.o cg_blas.o mmio.o 
OBJS=mmio.o cg.o util.o  second.o cg_blas.o 

all: conjugategradient

conjugategradient: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) $(CFLAGS) -o $@

#coverage study
coverage: 
	./conjugategradient matrix3x3.mtx vector3.mtx
	gcov *.c
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory out_html

clean:
	rm -Rf conjugategradient *.o *~ *.gcda *.gcov *.gcno out_html coverage.info gmon.out
