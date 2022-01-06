CC=mpicc
CFLAGS=-O3 -lm -Wall
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: gauss gauss_conv

gauss_conv: GaussSeidelSOR.c utils.c
	$(CC) $(CFLAGS) $(CONV) GaussSeidelSOR.c utils.c -o gauss_conv
gauss: GaussSeidelSOR.c utils.c
	$(CC) $(CFLAGS) GaussSeidelSOR.c utils.c -o gauss

clean:
	rm gauss gauss_conv


