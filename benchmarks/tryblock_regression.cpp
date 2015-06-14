#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <assert.h>
#include <unistd.h>

/*
#define TRYBLOCK      0
#define ALLREDUCE     0
#define IALLREDUCE    0
#define WITH_FAILURE  0
#define WITH_FAULT    0
*/

int main(int argc, char** argv)
{
	int rank, size;
	int i, j, k;
	double t1, t2, tf1 = 0.0, tf2 = 0.0;
	int* franks;
	int count;
	int rc;
	int numiters = 10000;
	int numskip  = 1000;
#ifdef WITH_FAILURE
	int numfailedskip = 1000;
#else
	int numfailedskip = 0;
#endif	/* WITH_FAILURE */
	int iamdead = false;

	MPI_Init(&argc, &argv);
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &size);

	MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN);
	
	MPI_Request tryreq, shrinkreq;
	MPI_Status trystat, shrinkstat;

	MPI_Timeout timeout;
	MPI_Timeout_set_seconds(&timeout, 1.0);

	if(0 == rank){
#ifdef TRYBLOCK
		printf("TRYBLOCK ");
#endif	/* TRYBLOCK */
		
#ifdef ALLREDUCE
		printf("ALLREDUCE ");
#endif	/* ALLREDUCE */

#ifdef IALLREDUCE
		printf("IALLREDUCE ");
#endif	/* IALLREDUCE */

#ifdef WITH_FAILURE
		printf("WITH_FAILURE ");
#endif	/* WITH_FAILURE */
		
#ifdef WITH_FAULT
		printf("WITH_FAULT ");
#endif	/* WITH_FAULT */
		printf("\n");
	}


#ifdef WITH_FAILURE
	for(i=1; i<argc; i++){
		if(rank == atoi(argv[i])){
			iamdead = true;
		}
	}
#endif	/* WITH_FAILRUE */


#ifdef WITH_FAULT
	int num_failed_div = 0;
	if(argc > 1)
		num_failed_div = atoi(argv[1]);
#endif // WITH_FAULT
	
//	if(rank == 0) printf("barrier done\n");
	
	for(i=0; i<numiters+numskip+numfailedskip; i++){

		if(i == numskip + numfailedskip)
			t1 = MPI_Wtime();


		int a[2], b[2];
#ifdef IALLREDUCE
		MPI_Request al_req;
		MPI_Iallreduce(a, b, 2, MPI_INT, MPI_BAND, world, &al_req);
		MPI_Wait(&al_req, MPI_STATUS_IGNORE);
#endif	/* IALLREDUCE */
		
#ifdef ALLREDUCE
		MPI_Allreduce(a, b, 2, MPI_INT, MPI_BAND, world);
#endif	/* ALLREDUCE */

#ifdef TRYBLOCK

#  ifdef WITH_FAILURE
		if(i == numskip){
			if(iamdead)
				*(int*)0 = 0;
			sleep(3);
			tf1 = MPI_Wtime();
		}
#  endif  /* WITH_FAILURE */
	
		MPI_Tryblock_start(world, MPI_TRYBLOCK_GLOBAL, &tryreq);

#  ifdef WITH_FAULT
		if(num_failed_div && !(rank % num_failed_div))
			MPI_Request_raise_error(tryreq, rank+size);
#  endif	/* WITH_FAULT */
		
		MPI_Tryblock_finish(tryreq, 0, NULL);
		MPI_Wait_local(&tryreq, &trystat, timeout);
		MPI_Request_free(&tryreq);

#  ifdef WITH_FAILURE
		if(i == numskip){
			tf2 = MPI_Wtime();
		}
#  endif  /* WITH_FAILURE */
		
#endif	/* TRYBLOCK */
		
	}	/* for */
	t2 = MPI_Wtime();

#ifdef WITH_FAILURE
	MPI_Comm newworld;

	MPI_Comm_ishrink(world, &newworld, &shrinkreq);
	assert(MPI_SUCCESS == MPI_Wait(&shrinkreq, &shrinkstat));

	MPI_Comm_free(&world);
	world = newworld;
	MPI_Comm_size(world, &size);
	MPI_Comm_rank(world, &rank);
#endif // WITH_FAILURE
	
	int print_rank = 0;

	double time = (t2-t1)/numiters;
	double max_time;
	MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_SUM, print_rank, world);

	double timef = tf2-tf1;
	double max_timef;
	MPI_Reduce(&timef, &max_timef, 1, MPI_DOUBLE, MPI_SUM, print_rank, world);

//	MPI_Allreduce(&time, &max_time, 1, MPI_DOUBLE, MPI_SUM, world);
	max_time /= size;			/* to get the average time not the maximum */
	max_timef /= size;			/* to get the average time not the maximum */

	if(rank == print_rank){
		printf("%d %.2f %.2f\n", size, max_time*1000000, max_timef*1000000);
		
	}
	MPI_Finalize();
	return 0;
}
