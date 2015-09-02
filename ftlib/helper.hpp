#ifndef _HELPER_HPP_
#define _HELPER_HPP_
#include <iostream>

#ifdef DEBUG
#define TRACE(r)												\
	{															\
		int rank;												\
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);					\
		if (rank == (r) || (r) < 0)								\
			printf("%d: %s:%d\n", rank, __func__, __LINE__);	\
	}
#else
#define TRACE(r)
#endif

namespace stencil{


}

#endif	/* _HELPER_HPP_ */
