BIN = ../
ODIR = obj
TARGET = $(BIN)/mpi-layer2
LIBS = -lm
CC = /usr/bin/mpicc
CFLAGS = -O2 -Wall

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, $(ODIR)/%.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

$(ODIR)/%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f $(ODIR)/*.o
	-rm -f $(TARGET)
