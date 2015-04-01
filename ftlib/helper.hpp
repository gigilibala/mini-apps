
#include <iostream>

/* public timing functions */
double gtime;
#define TICK()     (gtime = MPI_Wtime())
#define TOCK(time) (time += MPI_Wtime() - gtime)

namespace stencil{


}
