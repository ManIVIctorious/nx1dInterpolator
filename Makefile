
CC = gcc

CFLAGS += -w -Og -g

EXE = bin/nx1dInterpolator

OBJ += main.o
OBJ += nx1dInterpolator.o
OBJ += InputFunction.o

all: $(EXE) Makefile

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)
