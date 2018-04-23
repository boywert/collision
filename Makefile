CC=mpicc
CFLAGS= -O3 -I/users/home/boyd/local/include -L/users/home/boyd/local/lib  -lgsl -lgslcblas -lm

dsmc: hellomake.o hellofunc.o
     $(CC) -o dscm mpitest.c $(CFLAGS)
