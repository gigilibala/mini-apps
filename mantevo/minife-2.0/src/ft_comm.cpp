#include <ft_comm.hpp>
#include <iostream>
#include <ftlib.hpp>

namespace miniFE{


	FTComm* FTComm::ftcomm = NULL;

	void FTComm::init(bool spawned, int argc, char** argv){

		
		myargc = argc;
		myargv = argv;

		MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
		int rank, size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		MPI_Comm parent;
		MPI_Comm_get_parent(&parent);

		if(MPI_COMM_NULL == parent){
			/* this is the real world */
			spawned = false;
			restarted = false;
			MPI_Request dup_req;
			MPI_Status  dup_stat;
			/* duplicating communciator for now */
			MPI_Comm_idup(MPI_COMM_WORLD, &world, &dup_req);
			MPI_Wait(&dup_req, &dup_stat);

		}else{
			/* this is the spawned work */
			spawned = true;
			restarted = true;

			std::cout << "spawned, joining parent communicator" << std::endl;
			MPI_Intercomm_merge(parent, 1, &world);
			MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN);

			/* TODO: setup world */
//	    	MPI_Comm_free(&parent);
		}
	}


	void FTComm::repair(MPI_Request tb_req){
 
#if USING_FAMPI
		MPI_Comm shrinked_comm, spawned_comm, nwe;
		MPI_Comm local_world = world;
		
//		std::cout << "repairing the communicator" << std::endl;

		int err_code = MPI_ERR_PROC_FAILED;
		int spawn_size = 1;
/*
		MPI_Group fgroup;
		MPI_Get_failed_group(tb_req, 1, &err_code, local_world, &fgroup);
		MPI_Group_size(fgroup, &spawn_size);
		std::cout << "spawn size is:" << spawn_size << std::endl;
*/
		fampi_repair_comm_shrink(local_world, &shrinked_comm);
		MPI_Comm_free(&local_world);

//		std::cout << "going to spawn" << std::endl;
		fampi_repair_comm_spawn(shrinked_comm, spawn_size, myargc, myargv, &spawned_comm);
		MPI_Comm_free(&shrinked_comm);

		MPI_Intercomm_merge(spawned_comm, 1, &local_world);
		MPI_Comm_free(&spawned_comm);

		world = local_world;
		MPI_Comm_set_errhandler(world, MPI_ERRORS_RETURN);


		restarted = true;

		/* set size, argc, and argv stuff */
//		char** new_argv = new char*[myargc+1];
//		for(int i = 0; i < myargc; i++){
//			new_argv[i] = myargv[i];
//		}
//		new_argv[myargc-1] = (char*)"spawned=1";
//		new_argv[myargc]   = NULL;
//		fampi_repair_comm_spawn(shrinked_comm, 1, myargc+1, new_argv, &spawned_comm);
//		delete(new_argv);
#endif
	}

	FTComm* FTComm::get_instance(){
		if(NULL == ftcomm){
			ftcomm = new FTComm();
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

	bool FTComm::is_spawned(){
		return spawned;
	}

	bool FTComm::is_restarted(){
		return restarted;
	}

	void FTComm::set_restarted(){
		restarted = false;
	}

}
