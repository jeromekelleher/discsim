# Standard gcc/linux CFLAGS and LDFLAGS
CFLAGS=-g -Wall -Wextra
LDFLAGS=-lm -lgsl -lgslcblas
MAIN_LDFLAGS=${LDFLAGS} -lconfig 
ALL_SRC=sim.c util.c avl.c nystrom.c
CC=cc


main: main.c ${ALL_SRC} *.h
	${CC} ${CFLAGS} -o main main.c ${ALL_SRC} ${MAIN_LDFLAGS} 

tags:
	etags *.c *.h

clean:
	rm -f main
