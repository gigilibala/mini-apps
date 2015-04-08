#ifndef _Vector_hpp_
#define _Vector_hpp_

//@HEADER
// ************************************************************************
//
// MiniFE: Simple Finite Element Assembly and Solve
// Copyright (2006-2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// ************************************************************************
//@HEADER

#include <vector>
#include <string.h>
#include <assert.h>
#include <checkpoint.hpp>

namespace miniFE {


template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal>
struct Vector : ICheckpoint{
  typedef Scalar ScalarType;
  typedef LocalOrdinal LocalOrdinalType;
  typedef GlobalOrdinal GlobalOrdinalType;

  Vector(GlobalOrdinal startIdx, LocalOrdinal local_sz)
   : startIndex(startIdx),
     local_size(local_sz),
     coefs(local_size)
  {
    for(size_t i=0; i < local_size; ++i) {
    	coefs[i] = 0;
    }
  }

  ~Vector()
  {
  }

  GlobalOrdinal startIndex;
  LocalOrdinal local_size;
  std::vector<Scalar> coefs;

	int checkpoint(Checkpointer& cper){

		cper.cp_value<GlobalOrdinal>(startIndex);
		cper.cp_value<LocalOrdinal>(local_size);
		
		cper.cp_vector<Scalar>(coefs);

		return 0;
	}

	int restart(Checkpointer& cper){

		cper.r_value<GlobalOrdinal>(startIndex);
		cper.r_value<LocalOrdinal>(local_size);
		
		cper.r_vector<Scalar>(coefs);
		return 0;
	}
};

}//namespace miniFE

#endif

