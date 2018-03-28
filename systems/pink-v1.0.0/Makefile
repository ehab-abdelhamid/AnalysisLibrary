CC=mpicc
CFLAGS=-fopenmp
LIB=-lm -lrt -lgomp
PACKING_LIST=README \
Makefile \
pink.c \
shrink.c \
cluster.c \
cluster.h \
unionfind.c \
unionfind.h \
conv.c \
thresh.c \
verf.c

all:  pink shrink conv thresh

pink:  pink.o cluster.o unionfind.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

shrink:  shrink.o cluster.o unionfind.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

slink:  slink.o cluster.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

conv:  conv.o cluster.h
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

thresh:  thresh.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

verf:  verf.o cluster.o unionfind.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

slink.o:  clock.h

cluster.o:  cluster.h

unionfind.o:  unionfind.h

clean:
	rm -f pink shrink slink conv thresh verf *.o

tarball:
	tar -czvf pink.tar.gz $(PACKING_LIST)
