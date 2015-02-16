/* Copyright @ AMIN HASSANI 2015 */

#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpi-ext.h>

using namespace std;

typedef struct timeval timeval_t

/* public timing functions */
double gtime;
#define TICK()     (gtime = MPI_Wtime())
#define TOCK(time) (time += MPI_Wtime() - gtime)

int main(int argc, char *argv[])
{

	int rank, size;
	MPI_Comm world = MPI_COMM_WORLD;
	int matrix_size = 0;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &size);

	if(argc < 2){
		if(rank == 0)
			cout << "usage: stencil_1d.x <2d_matrix_size>" << endl;
		goto cleanup;
	}

	/* input processing */
	matrix_size = atoi(argv[1]);

	
	/* setup the random matrix */

	/* for loop for processing iterations*/

	/* recovery if needec */

cleanup:
	MPI_Finalize();
    return 0;
}


#if 0

//allocate memory from heap of size_row, size_col of char
char** memory_alloc(int size_row, int size_col){
	
  char** res = (char**)malloc(size_row*sizeof(char*));
  char* data = (char*) malloc(size_row*size_col*sizeof(char));

  if(res == NULL || data == NULL){
    printf("could not allocate memory\n");
    return NULL;
  }
  int i;
  for(i = 0; i < size_row; i++)
    res[i] = &data[i*size_col];
  return res;
}
/*prints the living matrix as output.
  zero shows cell is dead
  one shows cell is alive*/
void print_matrix(char** M, int size_row, int size_col, char mask){
  int i, j;
  for(i = 0; i < size_row; i++){
    for(j = 0; j < size_col; j++){
      printf("%d", M[i][j] & mask);
    }
    printf("\n");
  }
}

//initialize living matrix. border cells are all dead
void init_livings_matrix(char** M, int size_row, int size_col){
  //initializing border neighbours
  int i, j;
  for(i = 0; i < size_col; i++){
    M[0][i] = 0;
    M[size_row-1][i]=0;
  }
  for(i = 0; i < size_row; i++){
    M[i][0] = 0;
    M[i][size_col-1] = 0;
  }
	
  //initializing living cells
  for(i = 1; i <= size_row-2; i++){
    for(j = 1; j <= size_col-2; j++){
      M[i][j] = rand()%2;
    }
  }
}

//main function
//arguments: data_size max_generations num_threads
int main(int argc, char** argv){
  //printf("reached here\n");
  //check for input
  if(argc < 3){
    printf("not enough arguments\n");
    return -1;
  }
  int data_size, r_data_size;
  int size_row, size_col;
  int r_size_row, r_size_col; // real matrix size considering neighbours
  int max_generations, num_generations;
  int num_tasks, rank, rc;
  char **living_matrix;
  char **current_matrix;
  int changed, not_finished; // if 1 means next generations is changed so it needs another generation
  int i, j, k;
  int inc_size;
  char current_mask, next_mask;
  
  int *send_counts, *displacements;

  //MPI_Request* reqs;
  //MPI_Status* stats;

  MPI_Request requests[4];
  MPI_Status status[4];
  //simple arrays to find the place of neighbors
  static int neigs_i[] = {0, -1, -1, -1, 0, 1, 1 ,1};
  static int neigs_j[] = {1, 1, 0, -1, -1, -1, 0, 1};
  TIMESPEC t1, t2, t_res; // time take to execute

  rc = MPI_Init(&argc, &argv);
  if( rc != MPI_SUCCESS){
    printf("Initializing failed, Pleasae run program again\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Input argument processing
  data_size = atoi(argv[1]);
  max_generations = atoi(argv[2]);
  
  //Variable initialization
  r_data_size = data_size + 2;
  size_col = data_size;
  size_row = (int)ceil((float)data_size/num_tasks);
  if(rank == num_tasks-1 && data_size%size_row != 0) // the last processor gets less data than others
    size_row = data_size%size_row;
  r_size_row = size_row+2; 
  r_size_col = size_col+2;
  num_generations = 0;
  not_finished = 1;
  current_mask = 1;
  next_mask = 2;

  //Memory allocations
  //reqs = (MPI_Request*)malloc(num_tasks*sizeof(MPI_Request));
  //stats = (MPI_Status*)malloc(num_tasks*sizeof(MPI_Status));

  living_matrix = NULL;
  current_matrix = memory_alloc(r_size_row, r_size_col);
  send_counts = (int*)malloc(num_tasks*sizeof(int));
  displacements = (int*)malloc(num_tasks*sizeof(int));
  
  if(rank == 0){
    living_matrix = memory_alloc(r_data_size, r_data_size);
    srand(1);//can be changed to have different seeds.
    init_livings_matrix(living_matrix, r_data_size, r_data_size);

    //Initialize send buffer stuff
    for(i = 0; i < num_tasks; i++){
      send_counts[i] = r_size_row*r_size_col;
      displacements[i] = i*size_row*r_size_col;
    }
    int last_task_row;
    if((last_task_row = data_size%size_row) != 0){ // the last processor gets less data than others
      last_task_row += 2;
      send_counts[num_tasks-1] = r_size_col*last_task_row;
    }
  
  }else{
    living_matrix = memory_alloc(1,1);
  }

  //Scatter operations
  MPI_Scatterv(living_matrix[0], send_counts, displacements, MPI_CHAR, current_matrix[0], r_size_row*r_size_col, MPI_CHAR, 0, MPI_COMM_WORLD);
  
 //for gathering later
  for(i = 0; i < num_tasks; i++){
    send_counts[i] -= 2*r_size_col;
    displacements[i] += r_size_col;
    //printf("%d recv_counts: %d and displ: %d\n", i, send_counts[i], displacements[i]);
  }

  //Setting up the grid
  /*  int dims = num_tasks;
  int periods = 0;
  int reorder = 0;
  MPI_Comm cartcomm;
  MPI_Cart_create(MPI_COMM_WORLD, 1, &dims, &periods, reorder, &cartcomm);
  MPI_Comm_rank(cartcomm, &rank);
  MPI_Cart_coords(cartcomm, rank, 1, 
  */

  t1 = get_time();
  while(1){
    
    //num_iterations++;
    changed = 0;
    for(i = 1; i <= size_row; i++){
      for(j = 1; j <= size_col; j++){
	int alives = 0;
	for(k = 0; k < 8; k++){
	  int row = i + neigs_i[k];
	  int col = j + neigs_j[k];
	  if(current_matrix[row][col] & current_mask)
	    alives++;
	}
	//if current cell is alive
	if(current_matrix[i][j] & current_mask){
	  //will die of lonliness or starvation
	  if(alives <= 1 || alives >= 4){
	    current_matrix[i][j] &= ~next_mask;
	    changed = 1; // generation changed
	  }
	  else{//should stay alive
	    current_matrix[i][j] |= next_mask;
	  }
	}else{//if current cell is dead
	  if(alives == 3){
	    current_matrix[i][j] |= next_mask; // born the cell
	    changed = 1;
	  }
	  else
	    current_matrix[i][j] &= ~next_mask; // make sure its dead
	}
      }
    }
    //changing the bitmask order
    char tmp = current_mask;
    current_mask = next_mask;
    next_mask = tmp;
    
    num_generations++;
    not_finished = 0;
    MPI_Allreduce(&changed, &not_finished, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);

    if(!not_finished || num_generations == max_generations)
      break;

    //sending data to neighbours
    if(rank > 0)
      MPI_Isend(current_matrix[1], r_size_col, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, &requests[0]);
    if(rank < num_tasks-1)
      MPI_Isend(current_matrix[r_size_row-2], r_size_col, MPI_CHAR, rank+1, 1, MPI_COMM_WORLD, &requests[2]);
    if(rank > 0)
      MPI_Irecv(current_matrix[0], r_size_col, MPI_CHAR, rank-1, 1, MPI_COMM_WORLD, &requests[1]);
    if(rank < num_tasks-1)
      MPI_Irecv(current_matrix[r_size_row-1], r_size_col, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD, &requests[3]);
     
    if(num_tasks == 1) //in special case of one thread
      ;
    else if(rank == 0)
      MPI_Waitall(2, &requests[2], status);
    else if(rank == num_tasks-1)
      MPI_Waitall(2, &requests[0], status);
    else
      MPI_Waitall(4, &requests[0], status);

  }//end while
  //Gathering matrix back
  MPI_Gatherv(current_matrix[1], (r_size_row-2)*r_size_col, MPI_CHAR, living_matrix[0], send_counts, displacements, MPI_CHAR, 0, MPI_COMM_WORLD);  
  
  //do reduce to get correct time.
  //test if it needs barrier
  t2 = get_time();
  t_res = time_diff(t1, t2);
  int t_res_sec = t_res.tv_sec;
  
  MPI_Reduce(&(t_res.tv_sec), &t_res_sec, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  if(rank == 0){
    printf("Done with computation in %d secs and %d iters.\n", t_res_sec, num_generations);
  }
#if 0
  if(rank == 0){
    printf("final matrix:\n");
    print_matrix(living_matrix, r_data_size, r_data_size, (char)1);
  }
#endif

#if 0
  printf("my rank %d and current_generation:\n", rank);
  print_matrix(current_matrix, r_size_row, r_size_col, (char)1);
#endif

  MPI_Finalize();

  /*clean up and exit  */
  free(current_matrix[0]);
  free(current_matrix);
  free(living_matrix[0]);
  free(living_matrix);
  free(send_counts);
  free(displacements);

  exit(0);
  
  return 0;
}

#endif
