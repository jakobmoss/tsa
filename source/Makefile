# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis: Makefile for C program
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Settings
CC = gcc-5
CFLAGS = -Wall -std=gnu99
CFLAGS += -O3 -ffast-math -funroll-loops
# CFLAGS += -O2
CFLAGS += -fopenmp
LDLIBS += -lm -fopenmp

# Name of program and dependencies
NAME = powerspec
DEPEND = fileio.o arrlib.o tsfourier.o

# What to build
all: $(NAME)
	cp $(NAME) ../$(NAME).x

# Programs
$(NAME): $(NAME).o $(DEPEND)

# Housekeeping
clean:
	$(RM) $(NAME) *.o
