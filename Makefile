# Compiler
  CC   = gcc
# List of compiler flags
  CFLAGS += -g#                     # Enable debug symbols
  CFLAGS += -Og#                    # Set optimisation level, should be g if debug symbols are enabled
  CFLAGS += -march=native#          # Tune for current chipset, don't bother about backwards compatibility
 #CFLAGS += -mtune=native#          # Tune for current chipset, remain backwards compatible

  CFLAGS += -Werror#                # Treat warnings as errors
  CFLAGS += -Wall#                  # Enable base set of warnings
  CFLAGS += -Wextra#                # Enable additional warnings
  CFLAGS += -Wno-sign-compare#      # Disable sign-compare warning
  CFLAGS += -Wno-misleading-indentation#
  CFLAGS += -Wstrict-prototypes#    # Enable strict-prototypes warning  | not allowed in C++
  CFLAGS += -Wmissing-prototypes#   # Enable missing-prototypes warning | not allowed in C++
  CFLAGS += -Wno-unused-but-set-variable#
 #CFLAGS = -g -w#                   # Disable all warnings

# Resulting executable
  EXE = bin/nx1dInterpolator

# List of resulting object files
  OBJ += main.o
  OBJ += nx1dInterpolator.o
  OBJ += InputFunction.o

all: $(EXE)
# Build object files out of C-source files
%.o : %.c Makefile
	$(CC) $(CFLAGS) -c $<

# link all objects to create the executable
$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@


# allows to print out makefile variables, just type make print-VARIABLE
#  it will return VARIABLE = value_of_VARIABLE
print-%:
	@echo $* = $($*)

# remove all generated binary files
clean:
	rm -f $(OBJ) $(EXE)
