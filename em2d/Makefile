# GCC options
CC = gcc
CFLAGS = -std=c99 -pedantic -Wall -Wextra -Wno-stringop-truncation -O3

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic

# Clang options
#CC = clang
#CFLAGS = -Ofast

LDFLAGS = -lm


SOURCE = current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c

TARGET = zpic

all : clean $(TARGET)

$(TARGET) :
	$(CC) $(CFLAGS) $(SOURCE) $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET)
