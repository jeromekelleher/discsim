# Standard gcc/linux CFLAGS and LDFLAGS
CFLAGS=-g -Wall -Wextra
LDFLAGS=-lm -lgsl -lgslcblas
MAIN_LDFLAGS=${LDFLAGS} -lconfig 
ALL_SRC=sim.c 
#util.c kdtree.c torus.c
CC=cc


sim: main.c ${ALL_SRC} *.h
	${CC} ${CFLAGS} -o sim main.c ${ALL_SRC} ${MAIN_LDFLAGS} 

all_test: test_torus test_kdtree test_sim

check: all_test
	./test_torus
	./test_kdtree
	./test_sim


test_sim: test_sim.c ${ALL_SRC}
	${CC} ${CFLAGS} -o test_sim test_sim.c ${ALL_SRC} ${LDFLAGS} 

test_torus: test_torus.c ${ALL_SRC}
	${CC} ${CFLAGS} -o test_torus test_torus.c ${ALL_SRC} ${LDFLAGS} 

test_kdtree: test_kdtree.c ${ALL_SRC}
	${CC} ${CFLAGS} -o test_kdtree test_kdtree.c ${ALL_SRC} ${LDFLAGS} 

tags:
	etags *.c *.h

clean:
	rm -f sim test_torus test_kdtree test_sim
