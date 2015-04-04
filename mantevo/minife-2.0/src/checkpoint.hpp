#ifndef __checkpoint_hpp_
#define __checkpoint_hpp_


namespace miniFE{
	
class Checkpoint{

	Checkpoint(){/* do nothing */}
	~Checkpoint(){/* do nothing */}


	
	template<typename MatrixType, typename VectorType, typename Scalar>
	int checkpoint_write(int k, Scalar rtrans, Scalar oldrtrans, 
						 VectorType& p,
						 VectorType& r,
						 MatrixType& A,
						 char* buffer, int buf_size){

		int tmp_size;
		char* buf = buffer;
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
		return size;
	}
	

	template<typename MatrixType, typename VectorType, typename Scalar>
	int checkpoint_read(int& k, Scalar& rtrans, Scalar& oldrtrans, 
						VectorType& p,
						VectorType& r,
						MatrixType& A,
						char* buffer, int buf_size){
		
		int tmp_size;
		char* buf = buffer;

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

