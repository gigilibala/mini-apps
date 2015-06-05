# COPYRIGHT @ Amin Hassani 2015
CXX=mpicxx
CC=mpicc

#OBJS=$(patsubst %.cpp,%.o,$(SOURCES))
OBJS= main.o

CFLAGS   = -O2 -ggdb
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. -I../ftlib/
LDFLAGS  = 
FTLIB    = -L../ftlib -lfa_mpi
LIBS     = $(FTLIB)

USRPATH := $(shell echo $(USRPATH))

SOURCES=tryblock_regression.c

TARGETS=$(patsubst %.c,%.x,$(SOURCES))

all: $(TARGETS)

%.x: %.c
	$(CC) $(CFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

%.x: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

install: all
	cp $(TARGETS) $(USRPATH)/bin

clean: 
	rm $(TARGETS) *.o

.PHONY: all clean install 
