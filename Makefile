CC=mpicc
#CFLAGS= -O3 -I/users/home/boyd/local/include -L/users/home/boyd/local/lib  -lfftw3 -lgsl -lgslcblas -lm
CFLAGS= -O3 -lgsl -lgslcblas -lfftw3 -lm
all: fft dsmc dsmc_vhs
dsmc: mpitest.c
	$(CC) -o dsmc mpitest.c $(CFLAGS) 
dsmc_vhs: mpitest.c
	$(CC) -o dsmc_vhs mpitest.c $(CFLAGS) -DVHS
fft: collision.c
	$(CC) -o fft collision.c $(CFLAGS)
clean:
	rm -f dsmc dsmc_vhs fft
