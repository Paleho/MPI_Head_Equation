CC=mpicc
CFLAGS=-O3 -lm -Wall
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: redblack redblack_conv

redblack_conv: RedBlackSOR.c utils.c
	$(CC) $(CFLAGS) $(CONV) RedBlackSOR.c utils.c -o redblack_conv
redblack: RedBlackSOR.c utils.c
	$(CC) $(CFLAGS) RedBlackSOR.c utils.c -o redblack

clean:
	rm redblack redblack_conv


