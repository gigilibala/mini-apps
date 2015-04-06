#ifndef __checkpoint_hpp_
#define __checkpoint_hpp_

#include <fstream>
#include <vector>
#include <checkpoint.hpp>
#include <ft_comm.hpp>

#define CHECKPOINTING_TAG 123

namespace miniFE{

struct Checkpointer{

	
private:
	int remote_checkpoint_idx = 0;
	char* remote_checkpoint_buf[2];
	char* my_checkpoint_buf;
	int   my_checkpoint_size;
	int remote_rank;

	std::ostream* ofs;
	std::istream* ifs;
	
public:	
	void make_checkpoint_ready(){

		std::cout << "getting ready to checkpoint" << std::endl;
		MPI_Comm world = FTComm::get_instance()->get_world_comm();
		int rank;
		MPI_Comm_rank(world, &rank);
		char file_name[100];
		sprintf(file_name, "/tmp/%d.dat", rank);
		ofs = new std::ofstream(file_name, std::ofstream::binary);

	}

	void make_restart_ready(){

		std::cout << "getting ready to restart" << std::endl;
		MPI_Comm world = FTComm::get_instance()->get_world_comm();
		int rank;

		MPI_Comm_rank(world, &rank);
		char file_name[100];
		sprintf(file_name, "/tmp/%d.dat", rank);
		ifs = new std::ifstream(file_name, std::ofstream::binary);
	}

	Checkpointer(){


//		remote_checkpoint_idx = 0;
//		remote_checkpoint_buf[0] = (char*)malloc(max_size);
//		remote_checkpoint_buf[1] = (char*)malloc(max_size);
//		if(NULL == remote_checkpoint_buf[0] || NULL == remote_checkpoint_buf[1]){
//			MPI_Abort(FTComm::get_instance()->get_world_comm(), 0);
//		}
//		my_checkpoint_buf = (char*)malloc(max_size);
//		if(NULL == my_checkpoint_buf){
//			MPI_Abort(FTComm::get_instance()->get_world_comm(), 0);
//		}

		
	}
	~Checkpointer(){/* do nothing */}



	void checkpoint(){

		ofs->flush();
		((std::ofstream*)ofs)->close();
		
/*
		MPI_Comm world = FTComm::get_world_comm();
		MPI_Request reqs[2];
		MPI_Status  stats[2];
		int rank, size, remote;
		MPI_Comm_size(world, &size);
		MPI_Comm_rank(world, &rank);
		remote = (rank + (size/2)) % size;
		
		int recv_buf_idx = (remote_checkpoint_idx + 1) % 2;
		MPI_Irecv(remote_checkpoint_buf[recv_buf_idx], max_size, MPI_CHAR,
				  MPI_ANY_SOURCE, CHECKPOINTING_TAG, world, &reqs[0]);
		
		MPI_Isend(my_checkpoint_buf, my_checkpoint_size, MPI_CHAR, remote,
				  CHECKPOINTING_TAG, world, &reqs[1]);
		
		assert(MPI_SUCCESS == MPI_Waitall(2, reqs, stats));

		remote_rank = stats[0].MPI_SOURCE;
*/
	}

	void restart(){

		ifs->sync();
		((std::ifstream*)ifs)->close();
/*
		MPI_Comm world = FTComm::get_world_comm();
		MPI_Request reqs[2];
		MPI_Status  stats[2];
		int rank, size, remote;
		MPI_Comm_size(world, &size);
		MPI_Comm_rank(world, &rank);
		remote = (rank + (size/2)) % size;
		
		int recv_buf_idx = (remote_checkpoint_idx + 1) % 2;
		MPI_Irecv(remote_checkpoint_buf[recv_buf_idx], max_size, MPI_CHAR,
				  MPI_ANY_SOURCE, CHECKPOINTING_TAG, world, &reqs[0]);
		
		MPI_Isend(my_checkpoint_buf, my_checkpoint_size, MPI_CHAR, remote,
				  CHECKPOINTING_TAG, world, &reqs[1]);
		
		assert(MPI_SUCCESS == MPI_Waitall(2, reqs, stats));

		remote_rank = stats[0].MPI_SOURCE;
*/		
	}

	template<typename type>
	int cp_value(type& a){

		ofs->write((char*)&a, sizeof(type));
	}

	template<typename type>
	int cp_vector(std::vector<type>& v){
		int tmp_size = sizeof(type)*v.size();
		ofs->write((char*)&tmp_size, sizeof(int));
		ofs->write((char*)&v[0], tmp_size);
	}
	
	template<typename type>
	int r_value(type& a){
		ifs->read((char*)&a, sizeof(type));
	}

	template<typename type>
	int r_vector(std::vector<type>& v){
		int tmp_size = 0;
		ifs->read((char*)&tmp_size, sizeof(int));
		v.resize(tmp_size/sizeof(type));
		ifs->read((char*)&v[0], tmp_size);
	}


/*
	template<typename MatrixType, typename VectorType, typename Scalar>
	int checkpoint_write(int k, Scalar rtrans, Scalar oldrtrans, 
						 VectorType& p,
						 VectorType& r,
						 MatrixType& A){

		int tmp_size;
		char* buf = my_checkpoint_buf;
		int size = 0; 
		size += A.checkpoint_write(NULL, 0) +
			p.checkpoint_write(NULL, 0) +
			r.checkpoint_write(NULL, 0);
		
		size += 3 * sizeof(int);

		if(buf == NULL)
			return size;

		*(int*)buf = k; buf += sizeof(int);
		*(int*)buf = rtrans; buf += sizeof(int);
		*(int*)buf = oldrtrans; buf += sizeof(int);

		buf += A.checkpoint_write(buf, buf_size);
		buf += p.checkpoint_write(buf, buf_size);
		buf += r.checkpoint_write(buf, buf_size);

		assert(buf - buffer == size);

		my_checkpoint_size = size;
		
		checkpoint();
		
		return size;
	}
	

	template<typename MatrixType, typename VectorType, typename Scalar>
	int checkpoint_read(int& k, Scalar& rtrans, Scalar& oldrtrans, 
						VectorType& p,
						VectorType& r,
						MatrixType& A){
		
		int tmp_size;
		char* buf = my_checkpoint_buf;

		if(buf == NULL)
			return 0;

		k = *(int*)buf; buf += sizeof(int);
		rtrans = *(int*)buf; buf += sizeof(int);
		oldrtrans = *(int*)buf; buf += sizeof(int);

		buf += A.checkpoint_read(buf, buf_size);
		buf += p.checkpoint_read(buf, buf_size);
		buf += r.checkpoint_read(buf, buf_size);

		int size = checkpoint_write(k, rtrans, oldrtrans, p, r, A, buffer, buf_size);
		assert(buf - buffer == size);

		restart();
		
		return buf - buffer;
	}

*/
};
struct ICheckpoint{

public:
	virtual int checkpoint(Checkpointer& cper) = 0;
	virtual int restart(Checkpointer& cper) = 0;
};

}
#endif

