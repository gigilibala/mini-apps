#ifndef __ft_comm_hpp_
#define __ft_comm_hpp_

#include <mpi.h>
#include <mpi-ext.h>

namespace miniFE{

class FTComm{

private:
	static FTComm* ftcomm;
	MPI_Comm world;

	int myargc;
	char** myargv;

	bool spawned;
	bool restarted;

	/* call this right after the MPI_Init */
	FTComm(){spawned = false;}

	~FTComm(){}

public:

	void init(bool spawned, int argc, char** argv);

	void repair();

	static FTComm* get_instance();
	
	MPI_Comm get_world_comm();

	MPI_Comm shrink_and_spawn(MPI_Comm failed_comm);

	bool is_spawned();

	bool is_restarted();

	void set_restarted();

};
	
}

#endif
