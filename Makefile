CC=mpicc
#CFLAGS= -O3 -I/users/home/boyd/local/include -L/users/home/boyd/local/lib  -lgsl -lgslcblas -lm
CFLAGS= -O3 -lgsl -lgslcblas -lm
all: dsmc dsmc_vhs
dsmc: mpitest.c
	$(CC) -o dsmc mpitest.c $(CFLAGS) 
dsmc_vhs: mpitest.c
	$(CC) -o dsmc_vhs mpitest.c $(CFLAGS) -DVHS
clean:
	rm -f dsmc dsmc_vhs
