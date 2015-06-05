#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <assert.h>
#include <unistd.h>

/*
#define ALLREDUCE     0
#define IALLREDUCE    0
#define WITH_FAILURE  0
#define WITH_FAULT    0
*/

int main(int argc, char** argv)
{
	int rank, size;
	int i, j, k;
	double t1, t2;
	int* franks;
	int count;
	int rc;
	int numiters = 1000;
	int numskip  = 100;

	MPI_Init(&argc, &argv);
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &size);

	MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN);
	
	MPI_Request tryreq, shrinkreq;
	MPI_Status trystat, shrinkstat;

	MPI_Timeout timeout;
	MPI_Timeout_set_seconds(&timeout, 1.0);


#ifdef WITH_FAILURE
	for(i=1; i<argc; i++){
		if(rank == atoi(argv[i])){
			*(int*)0 = 0;
		}
	}

	sleep(1);
	MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tryreq);
retry:
	MPI_Tryblock_finish(tryreq, 0, NULL);
	rc = MPI_Wait_local(&tryreq, &trystat, timeout);
	if(MPI_SUCCESS != rc){
		if(MPI_ERR_TIMEOUT == rc){
			goto retry;
			
		}
		assert(MPI_ERR_TRYBLOCK_FOUND_ERRORS == rc);
	}
	MPI_Request_free(&tryreq);
#endif


#ifdef WITH_FAULT
	int num_failed_div = 0;
	if(argc > 1)
		num_failed_div = atoi(argv[1]);
#endif // WITH_FAULT
	
//	if(rank == 0) printf("barrier done\n");
	
	for(i=0; i<numiters+numskip; i++){
		if(i == numskip)
			t1 = MPI_Wtime();
//		printf("iter %d\n", i);

		int a[2], b[2];
#ifdef IALLREDUCE
		MPI_Request al_req;
		MPI_Iallreduce(a, b, 2, MPI_INT, MPI_BAND, world, &al_req);
		MPI_Wait(&al_req, MPI_STATUS_IGNORE);
#elif defined(ALLREDUCE)
		MPI_Allreduce(a, b, 2, MPI_INT, MPI_BAND, world);
#else
		MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tryreq);
#ifdef WITH_FAULT
		if(num_failed_div && !(rank % num_failed_div))
			MPI_Request_raise_error(tryreq, rank+size);
#endif
		
		MPI_Tryblock_finish(tryreq, 0, NULL);
		MPI_Wait_local(&tryreq, &trystat, timeout);
		MPI_Request_free(&tryreq);
#endif
	}		
	t2 = MPI_Wtime();

#ifdef WITH_FAILURE
	MPI_Comm newworld;
//	if(rank == 0) printf("calling shrink\n");
	MPI_Comm_ishrink(world, &newworld, &shrinkreq);
	assert(MPI_SUCCESS == MPI_Wait(&shrinkreq, &shrinkstat));
//	if(rank == 0) printf("called shrink\n");
	MPI_Comm_free(&world);
	world = newworld;
	MPI_Comm_size(world, &size);
	MPI_Comm_rank(world, &rank);
//	if(rank == 0) printf("shrinked properly\n");
#endif // WITH_FAILURE
	
	double time = (t2-t1)/numiters;
	double max_time;
//	MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_SUM, 0, world);
	MPI_Allreduce(&time, &max_time, 1, MPI_DOUBLE, MPI_SUM, world);
	max_time /= size;			/* to get the average time not the maximum */

	if(rank == 0){
//		printf("avg time (seconds)\n");
		printf("%d %.2f\n", size, max_time*1000000);
		
	}
	MPI_Finalize();
	return 0;
}
