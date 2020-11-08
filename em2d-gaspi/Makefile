# GCC options
CC = mpicc
CFLAGS = -Wall -Wextra -std=c99 -Wpedantic -Wno-stringop-truncation
CFLAGS += -O3

# Intel icc compiler
# CC = icc
# CFLAGS = -restrict -Ofast -std=c99 -pedantic -xCORE-AVX512 -mtune=skylake -xHost

# degub flags
# CFLAGS += -g
# CFLAGS += -ggdb
# CFLAGS += -finstrument-functions

# gprof
# CFLAGS += -pg

GPI2_HOME = /home/bruno/GPI-2-MPI
INCLUDE_DIR += $(GPI2_HOME)/include
LIBRARY_DIR += $(GPI2_HOME)/lib64

# LIB += GPI2-dbg
LIB += GPI2

LIB += m
LIB += pthread

LDFLAGS = -lm

LDFLAGS += $(addprefix -L,$(LIBRARY_DIR))
CFLAGS += $(addprefix -I,$(INCLUDE_DIR))

SOURCE = current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c gaspi_aux.c

TARGET = zpic

all : clean zpic

zpic :
	$(CC) $(CFLAGS) $(SOURCE) $(LDFLAGS) -o $@ $(addprefix -l, $(LIB))

clean:
	rm -f .gpi2.*
	rm -f $(TARGET)
