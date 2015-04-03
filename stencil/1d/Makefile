# COPYRIGHT @ Amin Hassani 2015
TARGET=stencil_1d.x
CXX=mpicxx
CC=mpicc

#OBJS=$(patsubst %.cpp,%.o,$(SOURCES))
OBJS= main.o

CFLAGS   = -O0 
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. -I../../ftlib/
LDFLAGS  = 
FTLIB    = -L ../../ftlib -lfa_mpi
LIBS     = $(FTLIB)

USRPATH := $(shell echo $(USRPATH))

OPTIONS  = 


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS)  $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

install: all
	cp $(TARGET) $(USRPATH)/bin

clean: 
	rm $(TARGET) *.o

.PHONY: all clean install 
