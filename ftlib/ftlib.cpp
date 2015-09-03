#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <assert.h>
#include <fstream>

#include "ftlib.hpp"
#include "helper.hpp"

using namespace std;

/* return error string */
static char error_string[MPI_MAX_ERROR_STRING];
static inline char* get_error_string(int rc){
	int strlen;
	MPI_Error_string(rc, error_string, &strlen);
	return error_string;
}

double gtime;

/**
 * spawns new processes and returnes the newly created communicator
 */
int fampi_repair_comm_spawn(MPI_Comm new_comm, int num_new_procs, int argc, char** argv, MPI_Comm* comm){
	TRACE(0);
	int exit_code = 0, rc;
	MPI_Comm intercomm;
	int error_codes;
	int rank, size;
	MPI_Comm_rank(new_comm, &rank);
		
	rc = MPI_Comm_spawn(argv[0], &argv[1], num_new_procs, MPI_INFO_NULL,
						0, new_comm, &intercomm, &error_codes);

	if(MPI_SUCCESS != rc){
		char error_str[200];
		int size = 200;
		MPI_Error_string(rc, error_str, &size);
		std::cout << "spawn failed in rank " << rank << " "
				  << "with error " << error_str << std::endl;
		exit_code = rc;
	}
	/* spawn section */
	*comm = intercomm;
	return exit_code;
}
 
/* shrinking a communicator by removing failed communicators  */
int fampi_repair_comm_shrink(MPI_Comm fcomm, MPI_Comm* comm){
	TRACE(0);
	int             rc, exit_code;
	MPI_Request     shrink_req, tb_req;
	MPI_Status      shrink_stat, tb_stat;
	MPI_Comm        newcomm;
	MPI_Timeout     shrink_timeout, tb_timeout;
	int             max_retries = 3;
	int             shrink_wait_max = 20;
	int             rank, size;

	MPI_Comm_rank(fcomm, &rank);
	MPI_Comm_size(fcomm, &size);
	

	/* shrink section */
	MPI_Timeout_set_seconds(&tb_timeout, 1.0);
	MPI_Timeout_set_seconds(&shrink_timeout, 2.0);

retry:
	
	MPI_Tryblock_start(fcomm, MPI_TRYBLOCK_GLOBAL, &tb_req);

	MPI_Comm_ishrink(fcomm, &newcomm, &shrink_req);

wait_shrink:
	rc = MPI_Wait_local(&shrink_req, &shrink_stat, shrink_timeout);
	if(MPI_ERR_TIMEOUT == rc && shrink_wait_max-- > 0){
		goto wait_shrink;
	}
	shrink_wait_max = 20;
	
retry_tb:
	assert(MPI_SUCCESS ==
		   MPI_Tryblock_finish_local(tb_req, 1, &shrink_req, shrink_timeout));
	
	rc = MPI_Wait_local(&tb_req, &tb_stat, tb_timeout);

	if(MPI_ERR_TIMEOUT == rc){
		if(rank == 0) cout << "shrink timed out, will retry" << endl;
		goto retry_tb;
	}

 	if(MPI_ERR_TRYBLOCK_FOUND_ERRORS == rc){
		cout << "shrink failed with error: " << get_error_string(rc) << endl;
		int test;
		MPI_Test(&shrink_req, &test, &shrink_stat);
		if(test)
			cout << "shrink had no problem" << endl;
		else
			std::cout << "shrink had problems" << std::endl;
		if(--max_retries > 0){
			cout << "going to retry" << endl;
			MPI_Request_free(&tb_req);
			MPI_Request_free(&shrink_req);
			MPI_Comm_free(&newcomm);
			goto retry;
		}else{
			cout << "will not be retried. will return with error\n" << endl;
			*comm = NULL;
		}
	}

	*comm = newcomm;
	MPI_Request_free(&tb_req);
	MPI_Request_free(&shrink_req);
	
	return exit_code;
}

/* TryBlock */
void TryBlockManager::push() {
	TRACE(0);
	tryblocks.push_back(NULL);
}

MPI_Request TryBlockManager::pop() {
	TRACE(0);
	TryBlock* tb = tryblocks.back();
	tryblocks.pop_back();
	if (NULL == tb) {
		std::cout << "Tryblock is null" << std::endl;
		return NULL;
	}
	return tb->tryblock_request;
}

int TryBlockManager::tryblock_start(MPI_Comm comm, int flag) {
	TRACE(0);
	tryblocks_be[tryblocks.size()-1]->start_timing();
	if (tryblocks.size() < 1) {
    	std:cout << __func__ << ": you should push first." << std::endl;
	}
	
	TryBlock* tb = 	tryblocks.back();
	if (NULL != tb)
		free(tb);

	tryblocks.back() = new TryBlock(comm);

	int rc = MPI_Tryblock_start(comm, flag,
								&tryblocks.back()->tryblock_request);
	
	tryblocks_be[tryblocks.size()-1]->end_timing();
	return rc;
}

int TryBlockManager::tryblock_finish(double reqs_timeout) {
	TRACE(0);
	tryblocks_be[tryblocks.size()-1]->start_timing();
	MPI_Timeout timeout;
	MPI_Timeout_set_seconds(&timeout, reqs_timeout);

	TryBlock* tb = tryblocks.back();

	int rc = MPI_Tryblock_finish_local(tb->tryblock_request,
									   tb->requests.size(),
									   tb->requests.data(),
									   timeout);

	tryblocks_be[tryblocks.size()-1]->end_timing();
	return rc;
}

int TryBlockManager::wait_for_tryblock_finish(double tb_timeout) {
	TRACE(0);
	tryblocks_be[tryblocks.size()-1]->start_timing();
	MPI_Timeout timeout;
	MPI_Timeout_set_seconds(&timeout, tb_timeout);
	int rc = MPI_Wait_local(
		&tryblocks.back()->tryblock_request, MPI_STATUS_IGNORE, timeout);

	tryblocks_be[tryblocks.size()-1]->end_timing();
	return rc;
}

int TryBlockManager::tryblock_finish_and_wait(double reqs_timeout, double tb_timout) {
	TRACE(0);
	tryblocks_be[tryblocks.size()-1]->start_timing();
	MPI_Timeout timeout;
	MPI_Timeout_set_seconds(&timeout, reqs_timeout);

	TryBlock* tb = tryblocks.back();
	int rc;
	do {
		MPI_Tryblock_finish_local(tb->tryblock_request,
								  tb->requests.size(),
								  tb->requests.data(),
								  timeout);
		rc = MPI_Wait_local(
			&tb->tryblock_request, MPI_STATUS_IGNORE, timeout);
	} while (rc == MPI_ERR_TIMEOUT);

	tryblocks_be[tryblocks.size()-1]->end_timing();
	return rc;
}

void TryBlockManager::add_requests(int count, MPI_Request* reqs) {
	TRACE(0);
	tryblocks.back()->requests.reserve(
		tryblocks.back()->requests.size() + count);

	for (int i=0; i<count; i++) {
		tryblocks.back()->requests.push_back(reqs[i]);
	}
}

void TryBlockManager::get_requests(int* count, MPI_Request** reqs) {
	TRACE(0);
	*count = tryblocks.back()->requests.size();
	*reqs = tryblocks.back()->requests.data();
}

void TryBlockManager::repair_comm(
	int argc, char** argv, MPI_Comm comm, MPI_Comm* out_comm) {
	TRACE(0);
	TryBlock* tb = tryblocks.back();

	int rc = MPI_Wait_local(&tb->tryblock_request, MPI_STATUS_IGNORE,
							MPI_TIMEOUT_ZERO);

	int array_of_error_codes[] = { MPI_ERR_PROC_FAILED };
	MPI_Group fgroup;
	int fsize;
	MPI_Get_failed_group(tb->tryblock_request, 1, array_of_error_codes,
						 comm, &fgroup);

	MPI_Group_size(fgroup, &fsize);


	if (timing_started) shrink_be->start_timing();
	MPI_Comm new_comm;
	fampi_repair_comm_shrink(comm, &new_comm);
	if (timing_started) shrink_be->end_timing();

	if (timing_started) spawn_be->start_timing();
	MPI_Comm spawned_comm;
	fampi_repair_comm_spawn(new_comm, fsize, argc, argv, &spawned_comm);
	if (timing_started) spawn_be->end_timing();
	
	if (type_ == 0) {			/* shrink */
		if (timing_started) merge_be->start_timing();
		MPI_Comm merged_comm;
		MPI_Intercomm_merge(spawned_comm, 1, &merged_comm);
		if (timing_started) merge_be->end_timing();

		*out_comm = merged_comm;
	} else if (type_ == 1) {		/* shrink and rebuild */
		if (timing_started) rebuild_be->start_timing();
		MPI_Comm rebuild_comm;
		MPI_Request rebuild_req;
		int* ranks = new int[fsize];
		for (int i = 0; i < fsize; i++)
			ranks[i] = i;
		MPI_Comm_irebuild(
			comm, 0, spawned_comm, fsize, ranks, 0, &rebuild_comm, &rebuild_req);
		MPI_Wait_local(&rebuild_req, MPI_STATUS_IGNORE, MPI_TIMEOUT_ZERO);
		if (timing_started) rebuild_be->end_timing();

		*out_comm = rebuild_comm;
	}

	MPI_Comm_free(&new_comm);
	MPI_Comm_free(&spawned_comm);
	MPI_Comm_set_errhandler(*out_comm, MPI_ERRORS_RETURN);

}

void TryBlockManager::start_timing() {
	timing_started = true;
}

void TryBlockManager::print_stats() {
	std::cout << tryblocks_be[0]->to_string() << std::endl;
	std::cout << tryblocks_be[1]->to_string() << std::endl;
	std::cout << shrink_be->to_string() << std::endl;
	std::cout << spawn_be->to_string() << std::endl;
	std::cout << merge_be->to_string() << std::endl;
}

void TryBlockManager::init() {
	TRACE(0);
	failed_cycle_ = 2;
	is_failed_ = false;
	repeat_ = 0;
	i_am_dead_ = false;
}

int TryBlockManager::repeat(MPI_Comm* world) {
	TRACE(0);
	repeat_++;
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	MPI_Comm parent;
	MPI_Comm_get_parent(&parent);
	if (parent != MPI_COMM_NULL) {	/* merged */
		if (type_ == 0) {	/* merge */
			MPI_Intercomm_merge(parent, 1, world);
		} else if (type_ == 1) { /* join */
			MPI_Request join_req;
			MPI_Comm_irebuild_join(parent, world, &join_req);
			MPI_Wait_local(&join_req, MPI_STATUS_IGNORE, MPI_TIMEOUT_ZERO);
		}

		is_failed_ = true;
		TRACE(0);
	} else {				/* from begining or repeat */
		if (repeat_ > 1) {	/* from repeat */
			TRACE(0);
			is_failed_ = true;
		} else {			/* from begining */
			/* It is from begining */
			/* Change the world communicator */
			TRACE(0);
			MPI_Comm_dup(MPI_COMM_WORLD, world);
			int rank;
			MPI_Comm_rank(*world, &rank);

			if (0 == rank) {
				process_file();
			}
			MPI_Bcast(failed_array_, fas_, MPI_INT, 0, *world);
		}
	}
	MPI_Comm_set_errhandler(*world, MPI_ERRORS_RETURN);
	return 0;
}

bool TryBlockManager::is_failed() {
	return is_failed_;
}

int TryBlockManager::failed_cycle() {
	return failed_cycle_;
}

bool TryBlockManager::am_i_dead(int cycle, int rank) {
	TRACE(0);
	int found_cycle = false;
	int found_rank  = false;
	int idx = 0;
	int cycles = failed_array_[idx++];
	for (int i = 0; i < cycles; i++) {
		int failed_cycle = failed_array_[idx++];
		int failed_ranks = failed_array_[idx++];
		if (failed_cycle == cycle)
			found_cycle = true;
		for (int j = 0; j < failed_ranks; j++) {
			if (rank == failed_array_[idx++])
				found_rank = true;
		}
	}
	return found_cycle && found_rank && !is_failed_;
}

/**
 * File structure for injecting failures.
 *
 * <num cycles failure>
 * <cycle1>  <num failed> <failed rank 1 > ... <failed rank (num failed)>
 * ...
 * <cyclen>  <num failed> <failed rank 1 > ... <failed rank (num failed)>
 */

int TryBlockManager::process_file() {
	std::ifstream ifs("ftlib.dat", std::ifstream::in);

	fas_ = 0;
	while (ifs.good()) {
		ifs >> failed_array_[fas_++];
	}
	ifs.close();

#if 1
	int i;
	for (i = 0; i < fas_; i++) {
		printf ("%d ", failed_array_[i]);
	}
	printf("\n");
#endif
}
