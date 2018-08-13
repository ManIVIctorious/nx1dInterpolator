# Compiler
  CC   = gcc
# List of compiler flags
  CFLAGS += -O2 -Wall -Wextra -Werror -march=native

# Resulting executable
  EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXE = $(EXEDIR)/nx1d-interpolator

# List of source files
  SRC += main.c
  SRC += nx1dInterpolator.c
  SRC += InputFunction.c

all: $(EXE)
# Build object files out of C-source files
$(SRC:.c=.o): $(SRC) Makefile
	$(CC) $(CFLAGS) $? -c

# link all objects to create the executable
$(EXE): $(SRC:.c=.o)
	$(CC) $(CFLAGS) $(LIB) $^ -o $@

# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(SRC:.c=.o) $(EXE)
