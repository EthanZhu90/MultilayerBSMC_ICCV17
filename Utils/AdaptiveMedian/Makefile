######################################################################
# Choose your favorite C compiler
CC = gcc

######################################################################
# -DNDEBUG prevents the assert() statements from being included in 
# the code.  If you are having problems running the code, you might 
# want to comment this line to see if an assert() statement fires.
FLAG1 = -DNDEBUG

######################################################################
# -DKLT_USE_QSORT forces the code to use the standard qsort() 
# routine.  Otherwise it will use a quicksort routine that takes
# advantage of our specific data structure to greatly reduce the
# running time on some machines.  Uncomment this line if for some
# reason you are unhappy with the special routine.
# FLAG2 = -DKLT_USE_QSORT

######################################################################
# Add your favorite C flags here.
CFLAGS = $(FLAG1) $(FLAG2)


LIB = -L/usr/local/lib -L/usr/lib

.SUFFIXES:  .c .o

all: adapmed.mexa64

.c.o:
	$(CC) -c $(CFLAGS) $<

adapmed.mexa64: adapmed.cpp 
	/usr/local/MATLAB/R2011b/bin/mex -I./ adapmed.cpp

clean:
	rm -f *.o 



