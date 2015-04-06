# Makefile for the fault-tolerant library used for FA-MPI sematics 
# COPYRIGHT @ AMIN HASSANI 2015

TARGET   = libfa_mpi
CXX      = mpicxx
CC       = mpicc

SOURCES  = $(shell find . -name '*.cpp' -print)

OBJS     = $(patsubst %.cpp,%.o,$(SOURCES))

CFLAGS   = -O0 -ggdb -fPIC
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
#	$(CXX) -shared -Wl,-soname,$(TARGET).so.1 -o $(TARGET).so.1.0 *.o
#	ld -G $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

install: all
	cp $(TARGET).a $(TARGET).so $(LIBDIR)
#	mv $(TARGET).so.1.0 $(LIBDIR)
#	ln -sf $(LIBDIR)/$(TARGET).so.1.0 $(LIBDIR)/$(TARGET).so.1
#	ln -sf $(LIBDIR)/$(TARGET).so.1   $(LIBDIR)/$(TARGET).so

clean: 
	rm $(TARGET).a $(TARGET).so *.o 

.PHONY: all clean install
