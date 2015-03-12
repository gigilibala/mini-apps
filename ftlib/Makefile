# Makefile for the fault-tolerant library used for FA-MPI sematics 
# COPYRIGHT @ AMIN HASSANI 2015

TARGET   = ftlib
CXX      = mpicxx
CC       = mpicc

SOURCES  = $(shell find . -name '*.cpp' -print)

OBJS     = $(patsubst %.cpp,%.o,$(SOURCES))

CFLAGS   = -O2 -g
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. 
LDFLAGS  =
LIBS     =

all: $(TARGET).a

$(TARGET).a: $(OBJS)
	$(AR) cr $@ $(OBJS)

$(TARGET).so: $(TARGET).a
	ld -G $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

install: all
	cp $(TARGET).a ../bin/

clean: 
	rm $(TARGET).a *.o

.PHONY: all clean install
