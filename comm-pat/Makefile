# Makefile for the fault-tolerant library used for FA-MPI sematics 
# COPYRIGHT @ AMIN HASSANI 2015

TARGET   = libprofiling
CXX      = mpicxx
CC       = mpicc

SOURCES  = $(shell find . -name '*.cpp' -print)

OBJS     = $(patsubst %.cpp,%.o,$(SOURCES))

CFLAGS   = -O2 -fPIC
# CFLAGS   = -O0 -ggdb -fPIC
# -Wall -Werror
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. 
LDFLAGS  =
LIBS     =

USRPATH := $(shell echo $(USRPATH))
LIBDIR   = $(USRPATH)/lib

# `echo $$USRPATH`/lib

all: $(TARGET).a $(TARGET).so

$(TARGET).a: $(OBJS)
	$(AR) -cvq $@ $(OBJS)
# cr

$(TARGET).so: $(OBJS)
	$(CXX) -shared -o $(TARGET).so *.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

install: all
	cp $(TARGET).a $(TARGET).so $(LIBDIR)

uninstall:
	rm $(LIBDIR)/$(TARGET).a $(LIBDIR)/$(TARGET).so 

clean: 
	rm $(TARGET).a $(TARGET).so *.o 

.PHONY: all clean install
