# Laptop options
NANOS6_HOME = /home/bruno/program-files/nanos6
INCLUDE_DIR += $(NANOS6_HOME)/include
LIBRARY_DIR += $(NANOS6_HOME)/lib64

CC=gcc
MCC=mcc
MPICC=mpicc
MPI_WRAPPERS=I_MPI_CC=$(MCC) MPICH_CC=$(MCC) OMPI_CC=$(MCC)

CFLAGS = -Wall -Wextra -std=c99 -Wpedantic -Wno-stringop-truncation
CFLAGS += -O3

# Intel icc compiler
# CC = icc
# MCC = imcc --ompss-2
# CFLAGS += --Wn,-restrict,-Ofast,-std=c99,-pedantic,-xCORE-AVX512,-mtune=skylake,-xHost

# Debug flags
# CFLAGS += -g
# CFLAGS += -ggdb
# CFLAGS += -finstrument-functions

# gprof
# CFLAGS += -pg

GPI2_HOME = ../GPI-2
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
	$(MCC) $(CFLAGS) $(SOURCE) $(LDFLAGS) -o $@ $(addprefix -l, $(LIB))

clean:
	rm -f .gpi2.*
	rm -f $(TARGET)
