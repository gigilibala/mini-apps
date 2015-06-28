# COPYRIGHT @ Amin Hassani 2015
CXX=mpicxx
CC=mpicc

CFLAGS   = -O2 -ggdb
CPPFLAGS = $(CFLAGS) 
CXXFLAGS = -I. -I../ftlib/
LDFLAGS  = -L../comm-pat
FTLIB    = -lprofiling
LIBS     = $(FTLIB)

DEFINES =
DEFINES += -DTRYBLOCK
# DEFINES += -DALLREDUCE
# DEFINES += -DIALLREDUCE
# DEFINES += -DWITH_FAILURE
# DEFINES += -DWITH_FAULT

USRPATH := $(shell echo $(USRPATH))

SOURCES=tryblock_regression.c ptp.c

TARGETS=$(patsubst %.c,%.x,$(SOURCES))

all: $(TARGETS)

%.x: %.c
	$(CC) $(CFLAGS) $(CXXFLAGS) $(DEFINES) $< -o $@ $(LDFLAGS) $(LIBS)

%.x: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEFINES) $< -o $@ $(LDFLAGS) $(LIBS)

install: all
	cp $(TARGETS) $(USRPATH)/bin

clean: 
	rm $(TARGETS)

.PHONY: all clean install 
