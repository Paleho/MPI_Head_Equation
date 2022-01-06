CC=mpicc
CFLAGS=-O3 -lm -Wall
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: jacobi_conv jacobi

jacobi_conv: Jacobi.c utils.c
	$(CC) $(CFLAGS) $(CONV) Jacobi.c utils.c -o jacobi_conv

jacobi: Jacobi.c utils.c
	$(CC) $(CFLAGS) Jacobi.c utils.c -o jacobi

clean:
	rm jacobi_conv jacobi
