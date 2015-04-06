#include <ft_comm.hpp>
#include <iostream>


namespace miniFE{


	FTComm* FTComm::ftcomm = NULL;

	void FTComm::init(bool spawned){

		if(NULL == ftcomm){
			ftcomm = new FTComm();
		
			if(!spawned){
				MPI_Request dup_req;
				MPI_Status  dup_stat;

				MPI_Comm_idup(FTComm::get_instance()->get_world_comm(), &ftcomm->world, &dup_req);
				MPI_Wait(&dup_req, &dup_stat);
			}else{

			}
		}
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
