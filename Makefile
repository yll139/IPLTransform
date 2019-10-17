TARGET = perf

OBJS = perf.o function.o matrix_operation.o
SOURCES = perf.c function.c matrix_operation.c

FLAG_I1 = -I ./include/mkl_include
FLAG_I2 = -I ./include/our_include

FLAG_L1 = -L ./include/mkl_lib_intel64
FLAG_l1 = -lmkl_rt

OPT 	= -O3

$(TARGET): $(OBJS)
	gcc -o perf $(OBJS) $(FLAG_L1) $(FLAG_l1)

perf.o: ./src/perf.c
	gcc -c ./src/perf.c -o perf.o $(FLAG_I1) $(FLAG_I2) $(OPT)

function.o: ./src/function.c
	gcc -c ./src/function.c -o function.o $(OPT) $(FLAG_I2)

matrix_operation.o: ./src/matrix_operation.c
	gcc -c ./src/matrix_operation.c -o matrix_operation.o $(OPT)

clean:
	rm -f *.o $(TARGET)
