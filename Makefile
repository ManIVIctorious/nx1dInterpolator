
CC = gcc

CFLAGS += -w -Og -g

EXE = bin/interpolator

OBJ += Interpolator.o
OBJ += InputFunction.o

all: $(EXE) Makefile

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@

clean:
	rm -f $(OBJ) $(EXE)
