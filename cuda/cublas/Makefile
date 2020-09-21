CC=nvcc
LD=${CC}
INCLUDES = -I$(OPENBLAS_ROOT)/include
CFLAGS+=-O3 $(INCLUDES) ${COVFLAGS}
DBGFLAGS = -pg -g

LDFLAGS+=-L$(OPENBLAS_ROOT)/lib -lopenblas 

#OBJS=blas.o cg.o util.o cg_blas.o mmio.o 
OBJS=mmio.o cg.o util.o  second.o cg_blas.o 

all: conjugategradient

conjugategradient: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -I$(OPENBLAS_ROOT)/include -o $@ -lcublas
	#$(LD) $(OBJS) $(LDFLAGS) $(DBGFLAGS) $(CFLAGS) -o $@

#coverage study
coverage: 
	./conjugategradient matrix3x3.mtx vector3.mtx
	gcov *.c
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory out_html

clean:
	rm -Rf conjugategradient *.o *~ *.gcda *.gcov *.gcno out_html coverage.info gmon.out

