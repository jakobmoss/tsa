# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis: Makefile for C program
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Settings
CC = gcc-5
CFLAGS = -Wall -std=gnu99
CFLAGS += -O3 -ffast-math -funroll-loops
CFLAGS += -fopenmp
LDLIBS += -lm -fopenmp

# Name of program and dependencies
NAME = powerspec
NAME2 = fclean
NAME3 = filter
DEPEND = fileio.o arrlib.o tsfourier.o window.o fmin.o pass.o

# What to build
all: $(NAME) $(NAME2) $(NAME3)
	cp $(NAME) ../$(NAME).x
	cp $(NAME2) ../$(NAME2).x
	cp $(NAME3) ../$(NAME3).x

# Programs
$(NAME): $(NAME).o $(DEPEND)

$(NAME2): $(NAME2).o $(DEPEND)

$(NAME3): $(NAME3).o $(DEPEND)


# Housekeeping
clean:
	$(RM) $(NAME) $(NAME2) $(NAME3) *.o
