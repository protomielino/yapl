CC      = gcc
CFLAGS  = -Wall -Wextra -O2
LDLIBS  = -lm -lraylib

SRC_COMMON = main.c maths.c sim.c

.PHONY: all grid kd run run-grid run-kd clean

all: grid kd

grid: $(SRC_COMMON) grid.c
	$(CC) $(CFLAGS) $^ -o yapl $(LDLIBS) -DUSE_GRID=1

kd: $(SRC_COMMON) kd.c
	$(CC) $(CFLAGS) $^ -o yapl_kd $(LDLIBS) -DUSE_GRID=0

run: run-grid

run-grid: grid
	./yapl

run-kd: kd
	./yapl_kd

clean:
	rm -f yapl yapl_kd
