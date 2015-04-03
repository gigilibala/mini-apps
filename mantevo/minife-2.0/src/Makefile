#-----------------------------------------------------------------------

MINIFE_TYPES =  \
        -DMINIFE_SCALAR=double   \
        -DMINIFE_LOCAL_ORDINAL=int      \
        -DMINIFE_GLOBAL_ORDINAL=int

MINIFE_MATRIX_TYPE = -DMINIFE_CSR_MATRIX
# MINIFE_MATRIX_TYPE = -DMINIFE_ELL_MATRIX

#-----------------------------------------------------------------------

CFLAGS = -O3 -g
CXXFLAGS = $(CFLAGS)

CPPFLAGS = -I. -I../utils ../fem $(MINIFE_TYPES) $(MINIFE_MATRIX_TYPE)
CPPFLAGS +=	-DHAVE_MPI
CPPFLAGS +=	-DUSING_FAMPI

#CPPFLAGS +=	-DMPICH_IGNORE_CXX_SEEK
#CPPFLAGS +=	-DMINIFE_REPORT_RUSAGE 

LDFLAGS=
LIBS=

CXX=mpicxx
CC=mpicc

#CXX=g++
#CC=gcc

include make_targets

