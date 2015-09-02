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
#include <helper.hpp>
#include <ftlib.hpp>

/* public timing functions */
extern double gtime;
#define TICK()     (gtime = MPI_Wtime())
#define TOCK(time) (time += MPI_Wtime() - gtime)

//#define DEBUG 1

#define WITH_FA 1

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
	MPI_Comm world, world_dup, new_world;
	int matrix_size = 0, rows, cols, internal_matrix_size, my_matrix_size;
	double time1, time2, time3;
	double *matrix = NULL, *tmp;
	double *neigh_top_ext, *neigh_bot_ext, *my_top_ext, *my_bot_ext;
	int iter = 0, max_iter;
	int top_neigh, bot_neigh;
	MPI_Request ee_reqs[4];
	MPI_Status  ee_stats[4];

	MPI_Request tb_req1,  tb_req2,  dup_req,  barrier_req,  shrink_req;
	MPI_Status  tb_stat1, tb_stat2, dup_stat, barrier_stat, shrink_stat;
	/* initiate the random number generator */
	srand(time(NULL));

	MPI_Timeout reqs_timeout, tb_timeout;
	bool iam_alive = true;
	bool reserved = false;

#if 0
	cout << "args:";
	for (int i = 0; i < argc; i++) {
		cout << " " << argv[i];
	}
	cout << endl;
#endif	

	if(argc < 3){
		if(rank == 0)
			cout << "usage: stencil_1d.x 0 <rows> <cols> <max_iter> <failed_rank1> <failed_rank2> ..." << endl;
		return 0;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	int spawned = atoi(argv[1]);
	if(!spawned){

		/* duplicating communciator for now */
		MPI_Comm_idup(MPI_COMM_WORLD, &world, &dup_req);
		MPI_Wait(&dup_req, &dup_stat);
		
	}else{
		MPI_Comm parent;
		MPI_Comm_get_parent(&parent);
		if(parent == MPI_COMM_NULL){
			std::cout << "no parent in spawned process about to abort" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 0);
		}
		MPI_Intercomm_merge(parent, 1, &world);
//		MPI_Comm_free(&parent);
		int recv_iter;
		MPI_Allreduce(&iter, &recv_iter, 1, MPI_INT, MPI_MAX,  world);
		iter = recv_iter + 1;
		MPI_Comm_rank(world, &rank);
		MPI_Comm_size(world, &size);

		char hostname[100];
		gethostname(hostname, 100);
		std::cout << "iter allreduced in child  " << iter << " "
				  << "in host: " << hostname << " "
				  << size << " " << rank << std::endl;

	}

//  fampi_repair_comm_spawn(world, size, argc, argv, &spawned_comm);
	
	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &size);
//	size-=3;

	/* input processing */
	rows = atoi(argv[2]);
	cols = atoi(argv[3]);

	max_iter = argc <= 4 ? 100 : atoi(argv[4]);
	for(int i=5; i<argc; i++){
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
	MPI_Timeout_set_seconds(&reqs_timeout, 1.0);
	
	/* make sure everyone is here */
//	MPI_Ibarrier(world, &barrier_req);
//	MPI_Wait(&barrier_req, &barrier_stat);

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
	for (; iter < max_iter; ++iter){

		if(0 == rank && 0 == iter%1)
			cout << "iter " << iter << ',' << rank << ']' << endl;
#if WITH_FA
		MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tb_req1);
#endif
		int req_i = 0;

		if(rank == 9 && (iter == 1 || iter == 3))
			*(int*)0 = 0;
		
		if(rank >= size) 
			goto skip;
		/* exchange externals */
		/* recv from top */
		MPI_Irecv(neigh_top_ext, cols, MPI_DOUBLE, top_neigh, 0, world, &ee_reqs[req_i++]);
		/* recv from bottom */
		MPI_Irecv(neigh_bot_ext, cols, MPI_DOUBLE, bot_neigh, 0, world, &ee_reqs[req_i++]);
		/* send to top */
		MPI_Isend(neigh_top_ext, cols, MPI_DOUBLE, top_neigh, 0, world, &ee_reqs[req_i++]);
		/* send to bottom */
		MPI_Isend(neigh_bot_ext, cols, MPI_DOUBLE, bot_neigh, 0, world, &ee_reqs[req_i++]);
		
	skip:
		/* compute */
		tmp = my_top_ext;
		for (int i = 0; i < my_matrix_size; ++i)
			heavy_op(tmp);

#if WITH_FA
	retry:
		assert(MPI_SUCCESS ==
			   MPI_Tryblock_finish_local(tb_req1, req_i, ee_reqs, reqs_timeout));
		rc = MPI_Wait_local(&tb_req1, &tb_stat1, tb_timeout);
		
		/* recovery if needed */
		if(MPI_ERR_TIMEOUT == rc){
			cout << rank << " tryblock timedout, have to retry" << endl;
			goto retry;
		}
		if(MPI_ERR_TRYBLOCK_FOUND_ERRORS == rc){

			int error_codes = MPI_ERR_PROC_FAILED;
			MPI_Comm comms[1];
			int comm_count;
			MPI_Get_failed_communicators(tb_req1, 1, &error_codes, 1, comms, &comm_count);

			if(comm_count > 0){
				if(rank == 0)
					cout << rank << " tryblock failed with " << comm_count <<
						" communicators failed, about to shrink" << endl;
				assert(comms[0] == world);
				
				MPI_Comm comm, spawned_comm, nwe;
				fampi_repair_comm_shrink(world, &comm);
				MPI_Comm_free(&world);

				argv[1][0] = '1'; /* shows it is spawned */
				
				fampi_repair_comm_spawn(comm, 1, argc, argv, &spawned_comm);
				MPI_Intercomm_merge(spawned_comm, 1, &world);
				MPI_Comm_free(&spawned_comm);
				
				argv[1][0] = '0';
				
				int recv_iter;
				MPI_Allreduce(&iter, &recv_iter, 1, MPI_INT, MPI_MAX,  world);
//				assert(iter == recv_iter);

				MPI_Comm_rank(world, &rank);
				MPI_Comm_size(world, &size);
				top_neigh = (rank-1+size)%size;
				bot_neigh = (rank+1)%size;
				if(rank == 0)
					std::cout << "iter allreduced in parent " << iter << " "
							  << size << " " << rank << std::endl;

			}
		}

		if(MPI_ERR_TIMEOUT == rc ||
		   MPI_ERR_TRANSIENT == rc ||
		   MPI_ERR_PERMANENT == rc)
			cout << "unhandled failure" << endl;
				
		for (int i = 0; i < req_i; ++i)
			MPI_Request_free(&ee_reqs[i]);
		MPI_Request_free(&tb_req1);
#else
		rc = MPI_Waitall(req_i, ee_reqs, MPI_STATUS_IGNORE);
#endif
	}
	TOCK(time1);
	MPI_Barrier(world);
//	if(0 == rank)
	cout << "total time is " << time1 << " in rank " << rank<< endl;
	
cleanup:

	if(matrix)
		free(matrix);
	sleep(2);					/* sleep for 2 seconds and then abort */
	MPI_Abort(world, 0);
	MPI_Finalize();
	return 0;
}
