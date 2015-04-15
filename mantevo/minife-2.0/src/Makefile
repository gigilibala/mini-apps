#-----------------------------------------------------------------------

MINIFE_TYPES =  \
        -DMINIFE_SCALAR=double   \
        -DMINIFE_LOCAL_ORDINAL=int64_t      \
        -DMINIFE_GLOBAL_ORDINAL=int64_t

MINIFE_MATRIX_TYPE = -DMINIFE_CSR_MATRIX
# MINIFE_MATRIX_TYPE = -DMINIFE_ELL_MATRIX

#-----------------------------------------------------------------------

CFLAGS = -O3 -g
CXXFLAGS = $(CFLAGS)

CPPFLAGS = -I. -I../utils -I../fem -I../../../ftlib $(MINIFE_TYPES) $(MINIFE_MATRIX_TYPE)
CPPFLAGS +=	-DHAVE_MPI
CPPFLAGS +=	-DUSING_FAMPI
CPPFLAGS +=	-DDEATH_ITER=2

#CPPFLAGS +=	-DCHECKPOINTING
#CPPFLAGS +=	-DNONBLOCKING
#CPPFLAGS +=	-DMPICH_IGNORE_CXX_SEEK
#CPPFLAGS +=	-DMINIFE_REPORT_RUSAGE 

LDFLAGS=
FTLIB= -L../../../ftlib -lfa_mpi
#FTLIB=
LIBS=$(FTLIB)

CXX=mpicxx
CC=mpicc

USRPATH := $(shell echo $(USRPATH))


#CXX=g++
#CC=gcc

include make_targets

