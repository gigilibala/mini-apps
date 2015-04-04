#ifndef __checkpoint_hpp_
#define __checkpoint_hpp_


namespace miniFE{
	
class Checkpoint{

private:
	char* remote_checkpoint_buf[2];
	char* my_checkpoint_buf;
	int remote_rank;
	
public:	
	Checkpoint(int max_size){


		int remote_checkpoint_idx = 0;
		remote_checkpoint_buf[0] = malloc(max_size);
		remote_checkpoint_buf[1] = malloc(max_size);
		if(NULL == remote_checkpoint_buf[0] || NULL == remote_checkpoint_buf[1]){
			MPI_Abort(MPI_COMM_WORLD, 0);
		}

		my_checkpoint_buf[1] = malloc(max_size);
		if(NULL == my_checkpoint_buf[0]){
			MPI_Abort(MPI_COMM_WORLD, 0);
		}
		
	}
	~Checkpoint(){/* do nothing */}

	
	void checkpoint(){

		MPI_Comm world = FTComm::get_world_comm();
		int rank, size, remote;
		MPI_Comm_size(world, &size);
		MPI_Comm_rank(world, &rank);

		remote = (rank + (size/2)) % size;

		MPI_Irecv();
		MPI_Isend();
		
	}

	void restart(){

	}
	
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

		send_checkpoint();
		
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
		return buf - buffer;
	}
};

}
#endif

