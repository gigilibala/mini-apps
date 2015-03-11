/* Copyright @ AMIN HASSANI 2015 */
/* make everything nonblocking */
#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <assert.h>
#include "../helper/helper.hpp"


//#define DEBUG 1

static inline void mult_div_op(double *value)
{
	(*value) *= (*value);
	if (*value) (*value) /= (*value);
}

#define heavy_op mult_div_op

using namespace std;

int main(int argc, char *argv[])
{
	int rc;
	int rank, size;
	MPI_Comm world = MPI_COMM_WORLD, world_dup, new_world;
	int matrix_size = 0, rows, cols, internal_matrix_size, my_matrix_size;
	double time1, time2, time3;
	double *matrix = NULL, *tmp;
	double *neigh_top_ext, *neigh_bot_ext, *my_top_ext, *my_bot_ext;
	int max_iter;
	int top_neigh, bot_neigh;
	MPI_Request ee_reqs[4];
	MPI_Status  ee_stats[4];

	MPI_Request tb_req1,  tb_req2,  dup_req,  barrier_req,  shrink_req;
	MPI_Status  tb_stat1, tb_stat2, dup_stat, barrier_stat, shrink_stat;
	/* initiate the random number generator */
	srand(time(NULL));

	MPI_Timeout reqs_timeout, tb_timeout;
	bool iam_alive = true;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &size);
	MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN);

	if(argc < 3){
		if(rank == 0)
			cout << "usage: stencil_1d.x <rows> <cols> <max_iter> <failed_rank1> <failed_rank2> ..." << endl;
		goto cleanup;
	}

	/* duplicating communciator for now */
	MPI_Comm_idup(world, &world_dup, &dup_req);
	MPI_Wait(&dup_req, &dup_stat);
	
	/* input processing */
	rows = atoi(argv[1]);
	cols = atoi(argv[2]);

	max_iter = argc <= 3 ? 100 : atoi(argv[3]);
	for(int i=4; i<argc; i++){
		if(rank == atoi(argv[i]))
			iam_alive = false;
	}
		
	
	/* create and setup the random matrix */
	matrix_size = (rows + 2) * cols;
	internal_matrix_size = (rows - 2) * cols;
	my_matrix_size = rows * cols;
	matrix = (double*)malloc(sizeof(double) * matrix_size);
	tmp = matrix;
	assert(tmp);
	assert(matrix);
	for (int i = 0; i < matrix_size; ++i, tmp++){
		*tmp = (double)rand();
	}
	neigh_top_ext = matrix;
	neigh_bot_ext = matrix + (rows + 1) * cols;
	my_top_ext    = matrix + cols;
	my_bot_ext    = matrix + rows * cols;
	

	/* setting up neighbors */
	top_neigh = (rank-1+size)%size;
	bot_neigh = (rank+1)%size;
	
	MPI_Timeout_set_seconds(&tb_timeout, 1.0);
	MPI_Timeout_set_seconds(&reqs_timeout, 2.0);
	
	/* make sure everyone is here */
	MPI_Ibarrier(world, &barrier_req);
	MPI_Wait(&barrier_req, &barrier_stat);

	if(0 == rank)
		cout << "barrier done!" << endl;

	if(!iam_alive){
		cout << "rank " << rank << " is going to die" << endl;
		*(int*)0 = 0;
	}
#ifdef DEBUG
	sleep(10);
#endif // DEBUG

	TICK();
	/* for loop for processing iterations*/
	for (int iter = 0; iter < max_iter; ++iter){

		if(0 == rank && 0 == iter%20)
			cout << "iter " << iter << endl;

		MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tb_req1);

		/* exchange externals */
		int req_i = 0;
		/* recv from top */
		MPI_Irecv(neigh_top_ext, cols, MPI_DOUBLE, top_neigh, 0, world, &ee_reqs[req_i++]);
		/* recv from bottom */
		MPI_Irecv(neigh_bot_ext, cols, MPI_DOUBLE, bot_neigh, 0, world, &ee_reqs[req_i++]);
		/* send to top */
		MPI_Isend(neigh_top_ext, cols, MPI_DOUBLE, top_neigh, 0, world, &ee_reqs[req_i++]);
		/* send to bottom */
		MPI_Isend(neigh_bot_ext, cols, MPI_DOUBLE, bot_neigh, 0, world, &ee_reqs[req_i++]);
		
		/* compute */
		tmp = my_top_ext;
		for (int i = 0; i < my_matrix_size; ++i)
			heavy_op(tmp);

		assert(MPI_SUCCESS ==
			   MPI_Tryblock_finish_local(tb_req1, 4, ee_reqs, reqs_timeout));
		rc = MPI_Wait_local(&tb_req1, &tb_stat1, tb_timeout);
		
		/* recovery if needed */
		if(MPI_SUCCESS != rc){

			switch(rc){
			case MPI_ERR_TRYBLOCK_FAILED:
				if(0 == rank)
					cout << "tryblock failed, about to shrink\n" << endl;
				
				MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tb_req2);
				MPI_Comm_ishrink(world, &new_world, &shrink_req);
				assert(MPI_SUCCESS ==
					   MPI_Tryblock_finish_local(tb_req2, 1, &shrink_req, reqs_timeout));
				rc = MPI_Wait_local(&tb_req2, &tb_stat2, tb_timeout);
				if(MPI_SUCCESS != rc){
					cout << "tryblock failed in shrink!" << endl;
					goto cleanup;
				}else{
					int rank1, size1;
					MPI_Comm_rank(new_world, &rank1);
					MPI_Comm_size(new_world, &size1);
					cout << "(" << rank << "," << rank1 << ") new tryblock with size " << size << endl;
				}
				break;
			case MPI_ERR_TIMEOUT:
				if(0 == rank)
					cout << "tryblock timedout\n" << endl;
				break;
			}
		} else {
			for (int i = 0; i < 4; ++i)
				MPI_Request_free(&ee_reqs[i]);
			MPI_Request_free(&tb_req1);
		}
	}
	
	TOCK(time1);
	if(0 == rank)
		cout << "total time is " << time1 << endl;
cleanup:

	if(matrix)
		free(matrix);
	MPI_Finalize();
	return 0;
}
