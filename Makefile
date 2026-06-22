CC      = gcc
CFLAGS  = -Wall -Wextra -O2
LDLIBS  = -lm -lraylib

SRC = main.c maths.c sim.c grid.c kd.c qt.c

.PHONY: all run clean run-kd run-qt

all: yapl

yapl: $(SRC)
	$(CC) $(CFLAGS) $^ -o yapl $(LDLIBS)

run: yapl
	./yapl

run-kd: yapl
	./yapl --backend kd

run-qt: yapl
	./yapl --backend qt

clean:
	rm -f yapl
