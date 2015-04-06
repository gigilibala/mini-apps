#include <ft_comm.hpp>
#include <iostream>
#include <ftlib.hpp>

namespace miniFE{


	FTComm* FTComm::ftcomm = NULL;

	void FTComm::init(bool spawned){

		if(NULL != ftcomm)
			return;
		ftcomm = new FTComm();
		MPI_Comm local_world = ftcomm->get_world_comm();
        //FTComm::get_instance()->get_world_comm();
		
		MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

		if(!spawned){

			MPI_Request dup_req;
			MPI_Status  dup_stat;
			/* duplicating communciator for now */
			MPI_Comm_idup(MPI_COMM_WORLD, &local_world, &dup_req);
			MPI_Wait(&dup_req, &dup_stat);
		
		}else{
			MPI_Comm parent;
  			MPI_Comm_get_parent(&parent);
			if(parent == MPI_COMM_NULL){
				std::cout << "no parent in spawned process about to abort" << std::endl;
				MPI_Abort(MPI_COMM_WORLD, 0);
			}
			MPI_Intercomm_merge(parent, 1, &local_world);

			/* TODO: setup world */
//	    	MPI_Comm_free(&parent);

		}
	}


	void FTComm::repair(){
 
		MPI_Comm shrinked_comm, spawned_comm, nwe;
		MPI_Comm local_world = ftcomm->get_world_comm();
		
		fampi_repair_comm_shrink(local_world, &shrinked_comm);

		MPI_Comm_free(&local_world);

		/* TODO: set size, argc, and argv stuff */
		fampi_repair_comm_spawn(shrinked_comm, 1, 0, NULL, &spawned_comm);

		MPI_Intercomm_merge(spawned_comm, 1, &local_world);

		MPI_Comm_free(&spawned_comm);

	}

	FTComm* FTComm::get_instance(){
		if(NULL == ftcomm){
			std::cout << "instance is null" << std::endl;
		}
		return ftcomm;
	}
	
	MPI_Comm FTComm::get_world_comm(){
		return world;
	}

	MPI_Comm FTComm::shrink_and_spawn(MPI_Comm failed_comm){


		/* at the end */
		return world;
	}
}
