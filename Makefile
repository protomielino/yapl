CC      = gcc
CFLAGS  = -Wall -Wextra -O2
LDLIBS  = -lm -lraylib

SRC = main.c maths.c sim.c grid.c kd.c

.PHONY: all run clean

all: yapl

yapl: $(SRC)
	$(CC) $(CFLAGS) $^ -o yapl $(LDLIBS)

run: yapl
	./yapl

run-kd: yapl
	./yapl --backend kd

clean:
	rm -f yapl
