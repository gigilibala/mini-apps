#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int a = 1;

	if(rank == 0) {
		MPI_Send(&a, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	} else {
		MPI_Recv(&a, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	MPI_Finalize();
    return 0;
}
