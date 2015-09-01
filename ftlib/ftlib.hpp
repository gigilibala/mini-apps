/*
 * Copyright (c) 2013-2014 University of Alabama at Birmingham
 *                         Department of Computer and Information Sciences
 *                         All rights reserved.
 * $COPYRIGHT$
 *
 * Author: Amin Hassani (ahassani@cis.uab.edu)
 *******************************************************************************
 * 
 * FT library to work with tryblocks
 */

#include <vector>
#include <mpi.h>

#include "benchmark.hpp"

/* Keep these functions for backward compatibility */
int fampi_repair_comm_shrink(MPI_Comm fcomm, MPI_Comm* comm);
int fampi_repair_comm_spawn(MPI_Comm new_comm, int vsize, int argc, char** argv, MPI_Comm* comm);


/* Example:
   MPI_Comm world = MPI_COMM_WORLD;
   TryBlockManager tb_manager();

   tb_manager.push();

   tb_manager.tryblock_start(flag);
   tb_manager.add_requests();
   
       tb_manager.push();
       tb_manager.tryblock_start(flag);
	   tb_manager.add_requests();
	   tb_manager.tryblock_finish();
	   tb_manager.wait_for_tryblock_finish();
	   req = tb_manager.pop();
	   

   tb_manager.add_requests(req);
   tb_manager.tryblock_finish();
   tb_manager.wait_for_tryblock_finish();
   req2 = tb_manager.pop();
   
 */


/* Container for a TryBlock. */
class TryBlock {
	friend class TryBlockManager;

public:
	TryBlock(MPI_Comm comm) : comm(comm) { };
	~TryBlock() { };

private:
	MPI_Comm comm;
	std::vector<MPI_Request> requests;
	MPI_Request tryblock_request;

};


/* Managing TryBlocks such as starting, finishing, and waiting for them to
 * finish. */
class TryBlockManager {

public:
	/* Push down to hierarchy */
	void push();

	MPI_Request pop();

	int tryblock_start(MPI_Comm comm, int flag);

	int tryblock_finish(double reqs_timeout);

	int wait_for_tryblock_finish(double tb_timeout);
	
	int tryblock_finish_and_wait(double reqs_timeout, double tb_timout);

	void add_requests(int count, MPI_Request* reqs);

	void get_requests(int* count, MPI_Request** reqs);

	void repair_comm(int argc, char** argv, MPI_Comm world, MPI_Comm* out_world);

	void start_timing();

	void print_stats();

	TryBlockManager() {
		timing_started = false;
		shrink_be = new BenchmarkEntry("shrink");
		spawn_be =  new BenchmarkEntry("spawn");
		merge_be =  new BenchmarkEntry("merge");
		tryblocks_be[0] = new BenchmarkEntry("tb_1");
		tryblocks_be[1] = new BenchmarkEntry("tb_2");
	};
	~TryBlockManager() { };

private:
	std::vector<TryBlock*> tryblocks;
	bool timing_started;

	BenchmarkEntry* shrink_be;
	BenchmarkEntry* spawn_be;
	BenchmarkEntry* merge_be;
	BenchmarkEntry* tryblocks_be[2];

};

