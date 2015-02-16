exec_file=stencil_1d.x
CXX=mpicxx
CC=mpicc

#OBJS=$(patsubst %.cpp,%.o,$(SOURCES))
OBJS= main.o

CFLAGS   = -O3 -g
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. 
LDFLAGS  =
LIBS     =


all: $(exec_file)

$(exec_file): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS)  $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

install: all
	cp $(exec_file) ../../bin/

clean: 
	rm $(exec_file) *.o

.PHONY: all clean install