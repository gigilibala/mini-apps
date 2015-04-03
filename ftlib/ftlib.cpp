#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <assert.h>

#include "ftlib.hpp"

using namespace std;

/* return error string */
static char error_string[MPI_MAX_ERROR_STRING];
static inline char* get_error_string(int rc){
	int strlen;
	MPI_Error_string(rc, error_string, &strlen);
	return error_string;
}


/*
 * spawns new processes and returnes the newly created communicator
 */

int fampi_repair_comm_spawn(MPI_Comm new_comm, int vsize, int argc, char** argv, MPI_Comm* comm){
	int exit_code = 0, rc;
	MPI_Comm intercomm;
	int error_codes;
	int rank, size;
	MPI_Comm_rank(new_comm, &rank);
	
	std::cout << "about to spawn " << std::endl;
	argv[1][0] = '1'; 			/* shows it is spawned */
	
	rc = MPI_Comm_spawn(argv[0], &argv[1], 1, MPI_INFO_NULL, 0, new_comm, &intercomm, &error_codes);

	if(MPI_SUCCESS != rc){
		char error_str[200];
		int size = 200;
		MPI_Error_string(rc, error_str, &size);
		std::cout << "spawn failed in rank " << rank << " "
				  << "with error " << error_str << std::endl;
	}
	/* spawn section */
	*comm = intercomm;
	return 0;
}
 
/* shrinking a communicator by removing failed communicators  */
int fampi_repair_comm_shrink(MPI_Comm fcomm, MPI_Comm* comm){

	int             rc;
	MPI_Request     shrink_req, tb_req;
	MPI_Status      shrink_stat, tb_stat;
	MPI_Comm        newcomm;
	MPI_Timeout     shrink_timeout, tb_timeout;
	int             max_retries = 3;

	int             rank, size;

	MPI_Comm_rank(fcomm, &rank);
	MPI_Comm_size(fcomm, &size);
	

	/* shrink section */
	MPI_Timeout_set_seconds(&tb_timeout, 1.0);
	MPI_Timeout_set_seconds(&shrink_timeout, 1.0);

retry:
	
	MPI_Tryblock_start(fcomm, MPI_TRYBLOCK_GLOBAL, &tb_req);

	MPI_Comm_ishrink(fcomm, &newcomm, &shrink_req);
	
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
	return 0;
}
