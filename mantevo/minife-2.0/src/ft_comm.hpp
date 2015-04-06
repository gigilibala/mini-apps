#ifndef __ft_comm_hpp_
#define __ft_comm_hpp_
#include <mpi.h>
#include <mpi-ext.h>

namespace miniFE{

class FTComm{

private:
	static FTComm* ftcomm;
	MPI_Comm world;

	/* call this right after the MPI_Init */
	FTComm(){

	}

	~FTComm(){

	}

public:

	static void init(bool spawned);

	static void repair();

	static FTComm* get_instance();
	
	MPI_Comm get_world_comm();

	MPI_Comm shrink_and_spawn(MPI_Comm failed_comm);
};
	
}

#endif
