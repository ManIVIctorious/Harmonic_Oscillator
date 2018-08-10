# Compiler and Compilerflags
  CC = gcc
# List of compiler flags
  CFLAGS = -O2 -Wall -Wextra -Werror -march=native

# Resulting executable
  EXEDIR = $(if ${MyLocalPath}, ${MyLocalPath}, bin)
  EXE = $(EXEDIR)/harmonic-oscillator


# List of linked libraries
  LIB = -lm `pkg-config --cflags --libs gsl`

# List of resulting object files
  OBJ += HarmonicOscillator.o

all: $(EXE)
# define rule to build object files out of C-source files
%.o: %.c
	$(CC) $(CFLAGS) -c $<

# link all objects to create the executable
$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIB) -o $@

clean:
	rm -f $(OBJ) $(EXE)
