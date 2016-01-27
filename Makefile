# Compiler
CC = gcc
# Compilerflags
CFLAGS = -g03 -Wall -Werror -Wstrict-prototypes -Wmissing-prototypes -mtune=native -I../include
#CFLAGS = -g03 -Wall -Wno-unused-but-set-variable -Werror -Wstrict-prototypes -Wmissing-prototypes -mtune=native -I../include
# Zu erstellende fertige Programme (mit relativem Pfad vorangestellt)
EXE = HarmonicOscillator.bin
# Objektdateien (.o bzw. .out), Librarys und Abhängigkeiten
# (wenn Abhängigkeit geändert -> make wird neu ausgeführt)
OBJ1 = main.o
LIB1 = -lm `pkg-config --cflags --libs gsl` #../include/differential.c
DEP1 = main.c Makefile #../include/header.h 
# Kommandoblock "all:" als Einsprungspunkt für make
all: $(EXE)
#
$(EXE): $(OBJ1) $(DEP1)
	$(CC) $(CFLAGS) $(OBJ1) $(LIB1) -o $@

clean:
	rm $(OBJ1) $(EXE)
