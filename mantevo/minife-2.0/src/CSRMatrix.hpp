#ifndef _CSRMatrix_hpp_
#define _CSRMatrix_hpp_

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

#include <cstddef>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <checkpoint.hpp>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace miniFE {

template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal>
struct
CSRMatrix : ICheckpoint{
  CSRMatrix()
   : has_local_indices(false),
     rows(), row_offsets(), row_offsets_external(),
     packed_cols(), packed_coefs(),
     num_cols(0)
#ifdef HAVE_MPI
     ,external_index(), external_local_index(), elements_to_send(),
      neighbors(), recv_length(), send_length(), send_buffer(), request()
#endif
  {
  }

  ~CSRMatrix()
  {}

  typedef Scalar        ScalarType;
  typedef LocalOrdinal  LocalOrdinalType;
  typedef GlobalOrdinal GlobalOrdinalType;

  bool                       has_local_indices;
  LocalOrdinal               num_cols;
  std::vector<GlobalOrdinal> rows;
  std::vector<LocalOrdinal>  row_offsets;
  std::vector<LocalOrdinal>  row_offsets_external;
  std::vector<GlobalOrdinal> packed_cols;
  std::vector<Scalar>        packed_coefs;

#ifdef HAVE_MPI
  std::vector<GlobalOrdinal> external_index;
  std::vector<GlobalOrdinal>  external_local_index;
  std::vector<GlobalOrdinal> elements_to_send;
  std::vector<int>           neighbors;
  std::vector<LocalOrdinal>  recv_length;
  std::vector<LocalOrdinal>  send_length;
  std::vector<Scalar>        send_buffer;
  std::vector<MPI_Request>   request;
#endif

	int checkpoint_size()
	{
		int size = 0;
		size += sizeof(bool);
		size += sizeof(LocalOrdinal);

		size += sizeof(GlobalOrdinal)* rows.size();
		size += sizeof(LocalOrdinal)*  row_offsets.size();
		size += sizeof(LocalOrdinal)*  row_offsets_external.size();
		size += sizeof(GlobalOrdinal)* packed_cols.size();
		size += sizeof(Scalar)*        packed_coefs.size();
		
#ifdef HAVE_MPI
		size += sizeof(GlobalOrdinal)* external_index.size();
		size += sizeof(GlobalOrdinal)* external_local_index.size();
		size += sizeof(GlobalOrdinal)* elements_to_send.size();
		size += sizeof(int)*           neighbors.size();
		size += sizeof(LocalOrdinal)*  recv_length.size();
		size += sizeof(LocalOrdinal)*  send_length.size();
//		size += sizeof(Scalar)*        send_buffer.size();
//		size += sizeof(MPI_Request)*   request.size();
#endif
		size += sizeof(int)*           11;
		return size;
	}
	
	int checkpoint(Checkpointer& cper){
		
		cper.cp_value<bool>(has_local_indices);
		cper.cp_value<LocalOrdinal>(num_cols);

		cper.cp_vector<GlobalOrdinal>(rows);
		cper.cp_vector<LocalOrdinal>(row_offsets);
		cper.cp_vector<LocalOrdinal>(row_offsets_external);
		cper.cp_vector<GlobalOrdinal>(packed_cols);
		cper.cp_vector<Scalar>(packed_coefs);
#ifdef HAVE_MPI
		cper.cp_vector<GlobalOrdinal>(external_index);
		cper.cp_vector<GlobalOrdinal>(external_local_index);
		cper.cp_vector<GlobalOrdinal>(elements_to_send);
		cper.cp_vector<int>(neighbors);
		cper.cp_vector<LocalOrdinal>(recv_length);
		cper.cp_vector<LocalOrdinal>(send_length);
#endif
		return 0;
	}

	int restart(Checkpointer& cper){

		cper.r_value<bool>(has_local_indices);
		cper.r_value<LocalOrdinal>(num_cols);

		cper.r_vector<GlobalOrdinal>(rows);
		cper.r_vector<LocalOrdinal>(row_offsets);
		cper.r_vector<LocalOrdinal>(row_offsets_external);
		cper.r_vector<GlobalOrdinal>(packed_cols);
		cper.r_vector<Scalar>(packed_coefs);
#ifdef HAVE_MPI
		cper.r_vector<GlobalOrdinal>(external_index);
		cper.r_vector<GlobalOrdinal>(external_local_index);
		cper.r_vector<GlobalOrdinal>(elements_to_send);
		cper.r_vector<int>(neighbors);
		cper.r_vector<LocalOrdinal>(recv_length);
		cper.r_vector<LocalOrdinal>(send_length);

#endif
		return 0;
	}

	void print_sizes(){
		std::cout << "has_local_indices: " << has_local_indices << std::endl;
		std::cout << "rows: " << rows.size() << std::endl;
		std::cout << "row_offsets: " << row_offsets.size() << std::endl;
		std::cout << "row_offsets_external: " << row_offsets_external.size() << std::endl;
		std::cout << "packed_cols: " << packed_cols.size() << std::endl;
		std::cout << "packed_coefs: " << packed_coefs.size() << std::endl;
		std::cout << "num_cols: " << num_cols << std::endl;
#ifdef HAVE_MPI
		std::cout << "external_index: " << external_index.size() << std::endl;
		std::cout << "external_local_index: " << external_local_index.size() << std::endl;
		std::cout << "elements_to_send: " << elements_to_send.size() << std::endl;
		std::cout << "neighbors: " << neighbors.size() << std::endl;
		std::cout << "recv_length: " << recv_length.size() << std::endl;
		std::cout << "send_length: " << send_length.size() << std::endl;
		std::cout << "send_buffer: " << send_buffer.size() << std::endl;
		std::cout << "request: " << request.size() << std::endl;
#endif
	}

  size_t num_nonzeros() const
  {
    return row_offsets[row_offsets.size()-1];
  }

  void reserve_space(unsigned nrows, unsigned ncols_per_row)
  {
    rows.resize(nrows);
    row_offsets.resize(nrows+1);
    packed_cols.reserve(nrows * ncols_per_row);
    packed_coefs.reserve(nrows * ncols_per_row);
  }

  void get_row_pointers(GlobalOrdinalType row, size_t& row_length,
                        GlobalOrdinalType*& cols,
                        ScalarType*& coefs)
  {
    ptrdiff_t local_row = -1;
    //first see if we can get the local-row index using fast direct lookup:
    if (rows.size() >= 1) {
      ptrdiff_t idx = row - rows[0];
      if (idx < rows.size() && rows[idx] == row) {
        local_row = idx;
      }
    }
 
    //if we didn't get the local-row index using direct lookup, try a
    //more expensive binary-search:
    if (local_row == -1) {
      typename std::vector<GlobalOrdinal>::iterator row_iter =
          std::lower_bound(rows.begin(), rows.end(), row);
  
      //if we still haven't found row, it's not local so jump out:
      if (row_iter == rows.end() || *row_iter != row) {
        row_length = 0;
        return;
      }
  
      local_row = row_iter - rows.begin();
    }

    LocalOrdinalType offset = row_offsets[local_row];
    row_length = row_offsets[local_row+1] - offset;
    cols = &packed_cols[offset];
    coefs = &packed_coefs[offset];
  }

};

}//namespace miniFE

#endif

