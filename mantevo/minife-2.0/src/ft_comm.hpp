#ifndef __ft_comm_hpp_
#define __ft_comm_hpp_

namespace miniFE{

class FTComm{

private:
	static MPI_Comm world;
	
public:
	/* call this right after the MPI_Init */
	FTComm(bool spawned){

		if(!spawned){
			MPI_Request dup_req;
			MPI_Status  dup_stat;

			MPI_Comm_idup(MPI_COMM_WORLD, &world, &dup_req);
			MPI_Wait(&dup_req, &dup_stat);

		}else{
			
		}
	}

	~FTComm(){

	}

	static MPI_Comm get_world_comm(){
		return world;
	}

	static MPI_Comm shrink_and_spawn(MPI_Comm& failed_comm){


		/* at the end */
		return world;
	}
		
};
	
}

#endif
