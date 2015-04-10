#ifndef _cg_solve_hpp_
#define _cg_solve_hpp_

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

#include <cmath>
#include <limits>
#include <unistd.h>

#include <Vector_functions.hpp>
#include <mytimer.hpp>
#include <checkpoint.hpp>

#include <outstream.hpp>

namespace miniFE {

int fault_num = 1;

template<typename Scalar>
void print_vec(const std::vector<Scalar>& vec, const std::string& name)
{
  for(size_t i=0; i<vec.size(); ++i) {
    std::cout << name << "["<<i<<"]: " << vec[i] << std::endl;
  }
}

template<typename VectorType>
bool breakdown(typename VectorType::ScalarType inner,
               const VectorType& v,
               const VectorType& w)
{
  typedef typename VectorType::ScalarType Scalar;
  typedef typename TypeTraits<Scalar>::magnitude_type magnitude;

//This is code that was copied from Aztec, and originally written
//by my hero, Ray Tuminaro.
//
//Assuming that inner = <v,w> (inner product of v and w),
//v and w are considered orthogonal if
//  |inner| < 100 * ||v||_2 * ||w||_2 * epsilon

  magnitude vnorm = std::sqrt(dot(v,v));
  magnitude wnorm = std::sqrt(dot(w,w));
  return std::abs(inner) <= 100*vnorm*wnorm*std::numeric_limits<magnitude>::epsilon();
}

template<typename OperatorType,
         typename VectorType,
         typename Matvec>
int
cg_solve(OperatorType& A,
         const VectorType& b,
         VectorType& x,
         Matvec matvec,
         typename OperatorType::LocalOrdinalType max_iter,
         typename TypeTraits<typename OperatorType::ScalarType>::magnitude_type& tolerance,
         typename OperatorType::LocalOrdinalType& num_iters,
         typename TypeTraits<typename OperatorType::ScalarType>::magnitude_type& normr,
         timer_type* my_cg_times)
{
  typedef typename OperatorType::ScalarType ScalarType;
  typedef typename OperatorType::GlobalOrdinalType GlobalOrdinalType;
  typedef typename OperatorType::LocalOrdinalType LocalOrdinalType;
  typedef typename TypeTraits<ScalarType>::magnitude_type magnitude_type;

  int exit_code = 0;
  timer_type t0 = 0, tWAXPY = 0, tDOT = 0, tMATVEC = 0, tMATVECDOT = 0;
  timer_type total_time = mytimer();
  timer_type tTRYBLOCK = 0, tCHECKPOINT = 0, tRECOVERY = 0;

  int myproc = 0;
  int mysize = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(FTComm::get_instance()->get_world_comm(), &myproc);
  MPI_Comm_size(FTComm::get_instance()->get_world_comm(), &mysize);
#ifdef USE_MPI_PCONTROL
  MPI_Pcontrol(1);
#endif
#endif

  if (!A.has_local_indices) {
    std::cerr << "miniFE::cg_solve ERROR, A.has_local_indices is false, needs to be true. This probably means "
       << "miniFE::make_local_matrix(A) was not called prior to calling miniFE::cg_solve."
       << std::endl;
    return exit_code;
  }

  size_t nrows = A.rows.size();
  LocalOrdinalType ncols = A.num_cols;

  VectorType r(b.startIndex, nrows);
  VectorType p(0, ncols);
  VectorType Ap(b.startIndex, nrows);

  normr = 0;
  magnitude_type rtrans = 0;
  magnitude_type oldrtrans = 0;

  LocalOrdinalType print_freq = max_iter/10;
  if (print_freq>50) print_freq = 50;
  if (print_freq<1)  print_freq = 1;
//  print_freq = 1;

  ScalarType one = 1.0;
  ScalarType zero = 0.0;

  TICK(); waxpby(one, x, zero, x, p); TOCK(tWAXPY);

//  print_vec(p.coefs, "p");

  TICK();
  matvec(A, p, Ap);
  TOCK(tMATVEC);

  TICK(); waxpby(one, b, -one, Ap, r); TOCK(tWAXPY);

  TICK(); rtrans = dot(r, r); TOCK(tDOT);


  normr = std::sqrt(rtrans);

  if (myproc == 0) {
    std::cout << "Initial Residual = "<< normr << std::endl;
  }

  magnitude_type brkdown_tol = std::numeric_limits<magnitude_type>::epsilon();

#ifdef MINIFE_DEBUG
  std::ostream& os = outstream();
  os << "brkdown_tol = " << brkdown_tol << std::endl;
#endif

  Checkpointer cper;
  MPI_Timeout reqs_timeout, tb_timeout;
  MPI_Timeout_set_seconds(&reqs_timeout, 0.1);
  MPI_Timeout_set_seconds(&tb_timeout, 2.0);
  int death_iter = 2;
  double tt1=0.0, tt2= 0.0;

  for(LocalOrdinalType k=1; k <= max_iter && normr > tolerance; ++k) {
restart:
	/* things we need to save here:
	   k, rtrans, oldrtrans, p, r, A
	*/
	/* checkpointing */
	  if(FTComm::get_instance()->is_restarted()){
	      TICK();
		  /* restarting from failure */
		  /* make checkpoint buffers or files ready */
#ifdef CHECKPOINTING
		  cper.make_restart_ready();
		  /* read the checkpoint */
		  cper.r_value(k);
		  cper.r_value(rtrans);
		  cper.r_value(oldrtrans);
		  p.restart(cper);
		  r.restart(cper);
		  A.restart(cper);
		  /* do the actual restarting */
		  cper.restart();
#else
		  k = death_iter + 1;
#endif
//		  k = 3;
		  FTComm::get_instance()->set_restarted();
		  TOCK(tCHECKPOINT);

	  }else if(0 == (k % cper.checkpoint_rate)){
		  /* make checkpoint buffers or files ready */
	      TICK();
#ifdef CHECKPOINTING
		  cper.make_checkpoint_ready();
		  /* write the checkpoint */
		  cper.cp_value(k);
		  cper.cp_value(rtrans);
		  cper.cp_value(oldrtrans);
		  p.checkpoint(cper);
		  r.checkpoint(cper);
		  A.checkpoint(cper);
		  /* do the actual checkpointing */
		  cper.checkpoint();
#endif
		  TOCK(tCHECKPOINT);
	  }

    if (k == 1) {
      TICK(); waxpby(one, r, zero, r, p); TOCK(tWAXPY);
    }
    else {
      oldrtrans = rtrans;
      TICK(); rtrans = dot(r, r); TOCK(tDOT);
      magnitude_type beta = rtrans/oldrtrans;
      TICK(); waxpby(one, r, beta, p, p); TOCK(tWAXPY);
    }

    normr = std::sqrt(rtrans);

    if (myproc == 0 && (k%print_freq==0 || k==max_iter)) {
      std::cout << "Iteration = "<<k<<"   Residual = "<<normr<<std::endl;
    }

    magnitude_type alpha = 0;
    magnitude_type p_ap_dot = 0;

#ifdef MINIFE_FUSED
    TICK();
    p_ap_dot = matvec_and_dot(A, p, Ap);
    TOCK(tMATVECDOT);
#else

	/* the place that magic exchage externals happen */
	/* Here put the tryblock start. does not need to go before this. Do the data
	 * exchange and checkpointing here too. */

#ifndef USING_FAMPI
    TICK(); matvec(A, p, Ap); TOCK(tMATVEC);
#else	
	MPI_Request tb_req;
	int done, rc, num_tb_retries = 3;

	MPI_Comm world = FTComm::get_instance()->get_world_comm();
	TICK(); MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tb_req); TOCK(tTRYBLOCK);

    TICK(); matvec(A, p, Ap); TOCK(tMATVEC);

	/* inject failure */
#if 1
	if(k == death_iter){
    	if(myproc == 2 /* || myproc == 10 */){
	        *(int*)0 = 0;
        }
    }
#endif

  retry_tryblocK:
	TICK();
//	MPI_Tryblock_finish_local(tb_req, A.request.size(), &A.request[0], reqs_timeout);
//	rc = MPI_Wait_local(&tb_req, MPI_STATUS_IGNORE, tb_timeout);

	MPI_Tryblock_finish(tb_req, A.request.size(), &A.request[0]);
//	rc = MPI_Wait_local(&tb_req, MPI_STATUS_IGNORE, MPI_TIMEOUT_ZERO);
	rc = MPI_Wait_local(&tb_req, MPI_STATUS_IGNORE, tb_timeout);

	TOCK(tTRYBLOCK);

	if(MPI_ERR_TIMEOUT == rc){
		if(num_tb_retries-- > 0)
			goto retry_tryblocK;
		else
			MPI_Abort(world, 0);
	}else if(MPI_SUCCESS != rc){
#if 1

		if(MPI_ERR_TRYBLOCK_FOUND_ERRORS == rc){

			TICK();
			int error_codes = MPI_ERR_PROC_FAILED;
			MPI_Comm comms[1];
			int comm_count;
			MPI_Get_failed_communicators(tb_req, 1, &error_codes, 1, comms, &comm_count);

			if(comm_count > 0){
				fault_num = 1;

/*
				if(rank == 0)
					cout << rank << " tryblock failed with " << comm_count <<
						" communicators failed, about to shrink" << endl;
*/
//				sleep(1);
				assert(comms[0] == FTComm::get_instance()->get_world_comm());
				FTComm::get_instance()->repair(tb_req);

				MPI_Request_free(&tb_req);
				for(int i = 0; i < A.request.size(); i++){
					MPI_Request_free(&A.request[i]);
				}
				TOCK(tRECOVERY);
				/* set the timeings */
//				std::cout << "recovered" << std::endl;

				/* free the memory now */
				exit_code = -1;
				goto cleanup;
			}
			TOCK(tRECOVERY);
		}

		/* free the request in exchange externals */
		MPI_Request_free(&tb_req);
		for(int i = 0; i < A.request.size(); i++){
			MPI_Request_free(&A.request[i]);
		}
#endif
	}
#endif
	
    TICK(); p_ap_dot = dot(Ap, p); TOCK(tDOT);
#endif

#ifdef MINIFE_DEBUG
    os << "iter " << k << ", p_ap_dot = " << p_ap_dot;
    os.flush();
#endif
    if (p_ap_dot < brkdown_tol) {
      if (p_ap_dot < 0 || breakdown(p_ap_dot, Ap, p)) {
        std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"<<std::endl;
#ifdef MINIFE_DEBUG
        os << "ERROR, numerical breakdown!"<<std::endl;
#endif
        //update the timers before jumping out.
		goto cleanup;
//        my_cg_times[WAXPY] += tWAXPY;
//        my_cg_times[DOT] += tDOT;
//        my_cg_times[MATVEC] += tMATVEC;
//        my_cg_times[TOTAL] += mytimer() - total_time;
//        return 0;
      } 
      else brkdown_tol = 0.1 * p_ap_dot;
    }
    alpha = rtrans/p_ap_dot;
#ifdef MINIFE_DEBUG
    os << ", rtrans = " << rtrans << ", alpha = " << alpha << std::endl;
#endif

#ifdef MINIFE_FUSED
    TICK();
    fused_waxpby(one, x, alpha, p, x, one, r, -alpha, Ap, r);
    TOCK(tWAXPY);
#else
    TICK(); waxpby(one, x, alpha, p, x);
            waxpby(one, r, -alpha, Ap, r); TOCK(tWAXPY);
#endif

    num_iters = k;
  }

#ifdef HAVE_MPI
#ifdef USE_MPI_PCONTROL
  MPI_Pcontrol(0);
#endif
#endif

cleanup:
  my_cg_times[WAXPY] += tWAXPY;
  my_cg_times[DOT] += tDOT;
  my_cg_times[MATVEC] += tMATVEC;
  my_cg_times[MATVECDOT] += tMATVECDOT;
  my_cg_times[CHECKPOINT] += tCHECKPOINT;
  my_cg_times[RECOVERY] += tRECOVERY;
  my_cg_times[TRYBLOCK] += tTRYBLOCK;
  my_cg_times[TOTAL] += mytimer() - total_time;

  return exit_code;
}

}//namespace miniFE

#endif

