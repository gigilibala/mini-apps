# COPYRIGHT @ Amin Hassani 2015
TARGET=stencil_1d.x
CXX=mpicxx
CC=mpicc

#OBJS=$(patsubst %.cpp,%.o,$(SOURCES))
OBJS= main.o

CFLAGS   = -O2 -g
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. 
LDFLAGS  = 
FTLIB    = ../../ftlib/ftlib.a
LIBS     = $(FTLIB)

# MAKEFILES = ../../ftlib

all: $(FTLIB) $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS)  $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

install: all
	cp $(TARGET) ../../bin/

clean: 
	rm $(TARGET) *.o

.PHONY: all clean install ftlib
