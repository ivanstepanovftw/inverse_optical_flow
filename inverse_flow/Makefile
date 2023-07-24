CC=gcc
C2=g++
CFLAGS=-Wall -Wextra -Wno-unused -pedantic -O4

backward_flow: backward_flow.cpp iio.o
	$(C2) $(CFLAGS) -o backward_flow backward_flow.cpp iio.o -lpng -ljpeg -ltiff
	
	
iio.o: iio.c
	$(CC) -std=c99 $(CFLAGS) -DNDEBUG -D_GNU_SOURCE -c iio.c

clean:
	rm -f backward_flow iio.o
