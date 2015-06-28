#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

bool file_opened = false;
std::ofstream pfile;
static int myrank;
char filename[256];

#define printstr()								\
	do{											\
		pfile << __func__ << std::endl;			\
	}while(0);

int MPI_Init(int *argc, char ***argv){
	int rc = PMPI_Init(argc, argv);
	PMPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	sprintf(filename, "%s-%d.txt", (*argv)[0], myrank);
	pfile.open(filename);
	printstr();
	return rc;
}
int MPI_Finalize(void){
	printstr();
	pfile.flush();
	pfile.close();
	return PMPI_Finalize();
}
int MPI_Initialized(int *flag){
	return PMPI_Initialized(flag);
}
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
	printstr();
	return PMPI_Init_thread(argc, argv, required, provided);
}
int MPI_Abort(MPI_Comm comm, int errorcode){
	printstr();
	return PMPI_Abort(comm, errorcode);
}
int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win){
	printstr();
	return PMPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, op, win);
}
int MPI_Add_error_class(int *errorclass){
	printstr();
	return PMPI_Add_error_class(errorclass);
}
int MPI_Add_error_code(int errorclass, int *errorcode){
	printstr();
	return PMPI_Add_error_code(errorclass, errorcode);
}
int MPI_Add_error_string(int errorcode, const char *string){
	printstr();
	return PMPI_Add_error_string(errorcode, string);
}
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
int MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
}
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
}
int MPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
}
int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr){
	printstr();
	return PMPI_Alloc_mem(size, info, baseptr);
}
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
	printstr();
	return PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}
int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
}
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
int MPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
}
int MPI_Alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
}
int MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
}
int MPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm){
	printstr();
	return PMPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
}
int MPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ialltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, request);
}
int MPI_Barrier(MPI_Comm comm){
	printstr();
	return PMPI_Barrier(comm);
}
int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ibarrier(comm, request);
}
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
	printstr();
	return PMPI_Bcast(buffer, count, datatype, root, comm);
}
int MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
	printstr();
	return PMPI_Bsend(buf, count, datatype, dest, tag, comm);
}
int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ibcast(buffer, count, datatype, root, comm, request);
}
int MPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Bsend_init(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Buffer_attach(void *buffer, int size){
	printstr();
	return PMPI_Buffer_attach(buffer, size);
}
int MPI_Buffer_detach(void *buffer, int *size){
	printstr();
	return PMPI_Buffer_detach(buffer, size);
}
int MPI_Cancel(MPI_Request *request){
	printstr();
	return PMPI_Cancel(request);
}
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]){
	printstr();
	return PMPI_Cart_coords(comm, rank, maxdims, coords);
}
int MPI_Cart_create(MPI_Comm old_comm, int ndims, const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart){
	printstr();
	return PMPI_Cart_create(old_comm, ndims, dims, periods, reorder, comm_cart);
}
int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]){
	printstr();
	return PMPI_Cart_get(comm, maxdims, dims, periods, coords);
}
int MPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank){
	printstr();
	return PMPI_Cart_map(comm, ndims, dims, periods, newrank);
}
int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank){
	printstr();
	return PMPI_Cart_rank(comm, coords, rank);
}
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest){
	printstr();
	return PMPI_Cart_shift(comm, direction, disp, rank_source, rank_dest);
}
int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *new_comm){
	printstr();
	return PMPI_Cart_sub(comm, remain_dims, new_comm);
}
int MPI_Cartdim_get(MPI_Comm comm, int *ndims){
	printstr();
	return PMPI_Cartdim_get(comm, ndims);
}
int MPI_Close_port(const char *port_name){
	printstr();
	return PMPI_Close_port(port_name);
}
int MPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_accept(port_name, info, root, comm, newcomm);
}
int MPI_Comm_c2f(MPI_Comm comm){
	printstr();
	return PMPI_Comm_c2f(comm);
}
int MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode){
	printstr();
	return PMPI_Comm_call_errhandler(comm, errorcode);
}
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result){
	printstr();
	return PMPI_Comm_compare(comm1, comm2, result);
}
int MPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_connect(port_name, info, root, comm, newcomm);
}
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *function, MPI_Errhandler *errhandler){
	printstr();
	return PMPI_Comm_create_errhandler(function, errhandler);
}
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn, MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval, void *extra_state){
	printstr();
	return PMPI_Comm_create_keyval(comm_copy_attr_fn, comm_delete_attr_fn, comm_keyval, extra_state);
}
int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_create_group(comm, group, tag, newcomm);
}
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_create(comm, group, newcomm);
}
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval){
	printstr();
	return PMPI_Comm_delete_attr(comm, comm_keyval);
}
int MPI_Comm_disconnect(MPI_Comm *comm){
	printstr();
	return PMPI_Comm_disconnect(comm);
}
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_dup(comm, newcomm);
}
int MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request){
	printstr();
	return PMPI_Comm_idup(comm, newcomm, request);
}
int MPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_dup_with_info(comm, info, newcomm);
}
int MPI_Comm_free_keyval(int *comm_keyval){
	printstr();
	return PMPI_Comm_free_keyval(comm_keyval);
}
int MPI_Comm_free(MPI_Comm *comm){
	printstr();
	return PMPI_Comm_free(comm);
}
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag){
	printstr();
	return PMPI_Comm_get_attr(comm, comm_keyval, attribute_val, flag);
}
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int nodes[], const int degrees[], const int targets[], const int weights[], MPI_Info info, int reorder, MPI_Comm * newcomm){
	printstr();
	return PMPI_Dist_graph_create(comm_old, n, nodes, degrees, targets, weights, info, reorder,  newcomm);
}
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[], const int sourceweights[], int outdegree, const int destinations[], const int destweights[], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph){
	printstr();
	return PMPI_Dist_graph_create_adjacent(comm_old, indegree, sources, sourceweights, outdegree, destinations, destweights, info, reorder, comm_dist_graph);
}
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[], int maxoutdegree, int destinations[], int destweights[]){
	printstr();
	return PMPI_Dist_graph_neighbors(comm, maxindegree, sources, sourceweights, maxoutdegree, destinations, destweights);
}
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *inneighbors, int *outneighbors, int *weighted){
	printstr();
	return PMPI_Dist_graph_neighbors_count(comm, inneighbors, outneighbors, weighted);
}
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *erhandler){
	printstr();
	return PMPI_Comm_get_errhandler(comm, erhandler);
}
int MPI_Comm_get_info(MPI_Comm comm, MPI_Info *info_used){
	printstr();
	return PMPI_Comm_get_info(comm, info_used);
}
int MPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen){
	printstr();
	return PMPI_Comm_get_name(comm, comm_name, resultlen);
}
int MPI_Comm_get_parent(MPI_Comm *parent){
	printstr();
	return PMPI_Comm_get_parent(parent);
}
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group){
	printstr();
	return PMPI_Comm_group(comm, group);
}
int MPI_Comm_join(int fd, MPI_Comm *intercomm){
	printstr();
	return PMPI_Comm_join(fd, intercomm);
}
int MPI_Comm_rank(MPI_Comm comm, int *rank){
	printstr();
	return PMPI_Comm_rank(comm, rank);
}
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group){
	printstr();
	return PMPI_Comm_remote_group(comm, group);
}
int MPI_Comm_remote_size(MPI_Comm comm, int *size){
	printstr();
	return PMPI_Comm_remote_size(comm, size);
}
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val){
	printstr();
	return PMPI_Comm_set_attr(comm, comm_keyval, attribute_val);
}
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler){
	printstr();
	return PMPI_Comm_set_errhandler(comm, errhandler);
}
int MPI_Comm_set_info(MPI_Comm comm, MPI_Info info){
	printstr();
	return PMPI_Comm_set_info(comm, info);
}
int MPI_Comm_set_name(MPI_Comm comm, const char *comm_name){
	printstr();
	return PMPI_Comm_set_name(comm, comm_name);
}
int MPI_Comm_size(MPI_Comm comm, int *size){
	printstr();
	return PMPI_Comm_size(comm, size);
}
int MPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]){
	printstr();
	return PMPI_Comm_spawn(command, argv, maxprocs, info, root, comm, intercomm, array_of_errcodes);
}
int MPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[], const int array_of_maxprocs[], const MPI_Info array_of_info[], int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]){
	printstr();
	return PMPI_Comm_spawn_multiple(count, array_of_commands, array_of_argv, array_of_maxprocs, array_of_info, root, comm, intercomm, array_of_errcodes);
}
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_split(comm, color, key, newcomm);
}
int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm){
	printstr();
	return PMPI_Comm_split_type(comm, split_type, key, info, newcomm);
}
int MPI_Comm_test_inter(MPI_Comm comm, int *flag){
	printstr();
	return PMPI_Comm_test_inter(comm, flag);
}
int MPI_Compare_and_swap(void *origin_addr, void *compare_addr, void *result_addr, MPI_Datatype datatype, int target_rank, MPI_Aint target_disp, MPI_Win win){
	printstr();
	return PMPI_Compare_and_swap(origin_addr, compare_addr, result_addr, datatype, target_rank, target_disp, win);
}
int MPI_Dims_create(int nnodes, int ndims, int dims[]){
	printstr();
	return PMPI_Dims_create(nnodes, ndims, dims);
}
int MPI_Errhandler_free(MPI_Errhandler *errhandler){
	printstr();
	return PMPI_Errhandler_free(errhandler);
}
int MPI_Error_class(int errorcode, int *errorclass){
	printstr();
	return PMPI_Error_class(errorcode, errorclass);
}
int MPI_Error_string(int errorcode, char *string, int *resultlen){
	printstr();
	return PMPI_Error_string(errorcode, string, resultlen);
}
int MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
	printstr();
	return PMPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
}
int MPI_Fetch_and_op(void *origin_addr, void *result_addr, MPI_Datatype datatype, int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win){
	printstr();
	return PMPI_Fetch_and_op(origin_addr, result_addr, datatype, target_rank, target_disp, op, win);
}
int MPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iexscan(sendbuf, recvbuf, count, datatype, op, comm, request);
}
int MPI_Finalized(int *flag){
	printstr();
	return PMPI_Finalized(flag);
}
int MPI_Free_mem(void *base){
	printstr();
	return PMPI_Free_mem(base);
}
int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
	printstr();
	return PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
}
int MPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
}
int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm){
	printstr();
	return PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
}
int MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
}
int MPI_Get_address(const void *location, MPI_Aint *address){
	printstr();
	return PMPI_Get_address(location, address);
}
int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count){
	printstr();
	return PMPI_Get_count(status, datatype, count);
}
int MPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count){
	printstr();
	return PMPI_Get_elements(status, datatype, count);
}
int MPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count){
	printstr();
	return PMPI_Get_elements_x(status, datatype, count);
}
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win){
	printstr();
	return PMPI_Get(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win);
}
int MPI_Get_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype, void *result_addr, int result_count, MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win){
	printstr();
	return PMPI_Get_accumulate(origin_addr, origin_count, origin_datatype, result_addr, result_count, result_datatype, target_rank, target_disp, target_count, target_datatype, op, win);
}
int MPI_Get_library_version(char *version, int *resultlen){
	printstr();
	return PMPI_Get_library_version(version, resultlen);
}
int MPI_Get_processor_name(char *name, int *resultlen){
	printstr();
	return PMPI_Get_processor_name(name, resultlen);
}
int MPI_Get_version(int *version, int *subversion){
	printstr();
	return PMPI_Get_version(version, subversion);
}
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int index[], const int edges[], int reorder, MPI_Comm *comm_graph){
	printstr();
	return PMPI_Graph_create(comm_old, nnodes, index, edges, reorder, comm_graph);
}
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int index[], int edges[]){
	printstr();
	return PMPI_Graph_get(comm, maxindex, maxedges, index, edges);
}
int MPI_Graph_map(MPI_Comm comm, int nnodes, const int index[], const int edges[], int *newrank){
	printstr();
	return PMPI_Graph_map(comm, nnodes, index, edges, newrank);
}
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors){
	printstr();
	return PMPI_Graph_neighbors_count(comm, rank, nneighbors);
}
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[]){
	printstr();
	return PMPI_Graph_neighbors(comm, rank, maxneighbors, neighbors);
}
int MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges){
	printstr();
	return PMPI_Graphdims_get(comm, nnodes, nedges);
}
int MPI_Grequest_complete(MPI_Request request){
	printstr();
	return PMPI_Grequest_complete(request);
}
int MPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn, MPI_Grequest_cancel_function *cancel_fn, void *extra_state, MPI_Request *request){
	printstr();
	return PMPI_Grequest_start(query_fn, free_fn, cancel_fn, extra_state, request);
}
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result){
	printstr();
	return PMPI_Group_compare(group1, group2, result);
}
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup){
	printstr();
	return PMPI_Group_difference(group1, group2, newgroup);
}
int MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup){
	printstr();
	return PMPI_Group_excl(group, n, ranks, newgroup);
}
int MPI_Group_free(MPI_Group *group){
	printstr();
	return PMPI_Group_free(group);
}
int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup){
	printstr();
	return PMPI_Group_incl(group, n, ranks, newgroup);
}
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup){
	printstr();
	return PMPI_Group_intersection(group1, group2, newgroup);
}
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup){
	printstr();
	return PMPI_Group_range_excl(group, n, ranges, newgroup);
}
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup){
	printstr();
	return PMPI_Group_range_incl(group, n, ranges, newgroup);
}
int MPI_Group_rank(MPI_Group group, int *rank){
	printstr();
	return PMPI_Group_rank(group, rank);
}
int MPI_Group_size(MPI_Group group, int *size){
	printstr();
	return PMPI_Group_size(group, size);
}
int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2, int ranks2[]){
	printstr();
	return PMPI_Group_translate_ranks(group1, n, ranks1, group2, ranks2);
}
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup){
	printstr();
	return PMPI_Group_union(group1, group2, newgroup);
}
int MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ibsend(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message, MPI_Status *status){
	printstr();
	return PMPI_Improbe(source, tag, comm, flag, message, status);
}
int MPI_Imrecv(void *buf, int count, MPI_Datatype type, MPI_Message *message, MPI_Request *request){
	printstr();
	return PMPI_Imrecv(buf, count, type, message, request);
}
int MPI_Info_create(MPI_Info *info){
	printstr();
	return PMPI_Info_create(info);
}
int MPI_Info_delete(MPI_Info info, const char *key){
	printstr();
	return PMPI_Info_delete(info, key);
}
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo){
	printstr();
	return PMPI_Info_dup(info, newinfo);
}
int MPI_Info_free(MPI_Info *info){
	printstr();
	return PMPI_Info_free(info);
}
int MPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag){
	printstr();
	return PMPI_Info_get(info, key, valuelen, value, flag);
}
int MPI_Info_get_nkeys(MPI_Info info, int *nkeys){
	printstr();
	return PMPI_Info_get_nkeys(info, nkeys);
}
int MPI_Info_get_nthkey(MPI_Info info, int n, char *key){
	printstr();
	return PMPI_Info_get_nthkey(info, n, key);
}
int MPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag){
	printstr();
	return PMPI_Info_get_valuelen(info, key, valuelen, flag);
}
int MPI_Info_set(MPI_Info info, const char *key, const char *value){
	printstr();
	return PMPI_Info_set(info, key, value);
}
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm bridge_comm, int remote_leader, int tag, MPI_Comm *newintercomm){
	printstr();
	return PMPI_Intercomm_create(local_comm, local_leader, bridge_comm, remote_leader, tag, newintercomm);
}
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintercomm){
	printstr();
	return PMPI_Intercomm_merge(intercomm, high, newintercomm);
}
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status){
	printstr();
	return PMPI_Iprobe(source, tag, comm, flag, status);
}
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
}
int MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Irsend(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Issend(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Is_thread_main(int *flag){
	printstr();
	return PMPI_Is_thread_main(flag);
}
int MPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name){
	printstr();
	return PMPI_Lookup_name(service_name, info, port_name);
}
int MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status){
	printstr();
	return PMPI_Mprobe(source, tag, comm, message, status);
}
int MPI_Mrecv(void *buf, int count, MPI_Datatype type, MPI_Message *message, MPI_Status *status){
	printstr();
	return PMPI_Mrecv(buf, count, type, message, status);
}
int MPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Neighbor_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ineighbor_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
}
int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Neighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
}
int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ineighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
}
int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Neighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
}
int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ineighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
}
int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],  MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm){
	printstr();
	return PMPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls,  sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
}
int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
}
int MPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm){
	printstr();
	return PMPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
}
int MPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, request);
}
int MPI_Op_commutative(MPI_Op op, int *commute){
	printstr();
	return PMPI_Op_commutative(op, commute);
}
int MPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op){
	printstr();
	return PMPI_Op_create(function, commute, op);
}
int MPI_Open_port(MPI_Info info, char *port_name){
	printstr();
	return PMPI_Open_port(info, port_name);
}
int MPI_Op_free(MPI_Op *op){
	printstr();
	return PMPI_Op_free(op);
}
int MPI_Pack_external(const char datarep[], const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, MPI_Aint outsize, MPI_Aint *position){
	printstr();
	return PMPI_Pack_external(datarep, inbuf, incount, datatype, outbuf, outsize, position);
}
int MPI_Pack_external_size(const char datarep[], int incount, MPI_Datatype datatype, MPI_Aint *size){
	printstr();
	return PMPI_Pack_external_size(datarep, incount, datatype, size);
}
int MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm){
	printstr();
	return PMPI_Pack(inbuf, incount, datatype, outbuf, outsize, position, comm);
}
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size){
	printstr();
	return PMPI_Pack_size(incount, datatype, comm, size);
}
/*
int MPI_Pcontrol(const int level, ...){
	printstr();
	return PMPI_Pcontrol(level, ...);
}
*/
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status){
	printstr();
	return PMPI_Probe(source, tag, comm, status);
}
int MPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name){
	printstr();
	return PMPI_Publish_name(service_name, info, port_name);
}
int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win){
	printstr();
	return PMPI_Put(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win);
}
int MPI_Query_thread(int *provided){
	printstr();
	return PMPI_Query_thread(provided);
}
int MPI_Raccumulate(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request){
	printstr();
	return PMPI_Raccumulate(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, op, win, request);
}
int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Recv_init(buf, count, datatype, source, tag, comm, request);
}
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status){
	printstr();
	return PMPI_Recv(buf, count, datatype, source, tag, comm, status);
}
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
	printstr();
	return PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}
int MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
}
int MPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op){
	printstr();
	return PMPI_Reduce_local(inbuf, inoutbuf, count, datatype, op);
}
int MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
	printstr();
	return PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
}
int MPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
}
int MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
	printstr();
	return PMPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm);
}
int MPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ireduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm, request);
}
/*
int MPI_Register_datarep(const char *datarep, MPI_Datarep_conversion_function *read_conversion_fn, MPI_Datarep_conversion_function *write_conversion_fn, MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state){
	printstr();
	return PMPI_Register_datarep(datarep, read_conversion_fn, write_conversion_fn, dtype_file_extent_fn, extra_state);
}
*/
int MPI_Request_free(MPI_Request *request){
	printstr();
	return PMPI_Request_free(request);
}
int MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status){
	printstr();
	return PMPI_Request_get_status(request, flag, status);
}
int MPI_Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request){
	printstr();
	return PMPI_Rget(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_count, target_datatype, win, request);
}
int MPI_Rget_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype, void *result_addr, int result_count, MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request){
	printstr();
	return PMPI_Rget_accumulate(origin_addr, origin_count, origin_datatype, result_addr, result_count, result_datatype, target_rank, target_disp, target_count, target_datatype, op, win, request);
}
int MPI_Rput(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp, int target_cout, MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request){
	printstr();
	return PMPI_Rput(origin_addr, origin_count, origin_datatype, target_rank, target_disp, target_cout, target_datatype, win, request);
}
int MPI_Rsend(const void *ibuf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
	printstr();
	return PMPI_Rsend(ibuf, count, datatype, dest, tag, comm);
}
int MPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Rsend_init(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
	printstr();
	return PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
}
int MPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iscan(sendbuf, recvbuf, count, datatype, op, comm, request);
}
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
	printstr();
	return PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
}
int MPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
}
int MPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
	printstr();
	return PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
}
int MPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
}
int MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Send_init(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
	printstr();
	return PMPI_Send(buf, count, datatype, dest, tag, comm);
}
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,  MPI_Status *status){
	printstr();
	return PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm,  status);
}
int MPI_Sendrecv_replace(void * buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status *status){
	printstr();
	return PMPI_Sendrecv_replace( buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
}
int MPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
	printstr();
	return PMPI_Ssend_init(buf, count, datatype, dest, tag, comm, request);
}
int MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
	printstr();
	return PMPI_Ssend(buf, count, datatype, dest, tag, comm);
}
int MPI_Start(MPI_Request *request){
	printstr();
	return PMPI_Start(request);
}
int MPI_Startall(int count, MPI_Request array_of_requests[]){
	printstr();
	return PMPI_Startall(count, array_of_requests);
}
/*
int MPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status){
	printstr();
	return PMPI_Status_c2f(c_status, MPI_Ff_status);
}
int MPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status){
	printstr();
	return PMPI_Status_f2c(MPI_Ff_status, c_status);
}
*/
int MPI_Status_set_cancelled(MPI_Status *status, int flag){
	printstr();
	return PMPI_Status_set_cancelled(status, flag);
}
int MPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count){
	printstr();
	return PMPI_Status_set_elements(status, datatype, count);
}
int MPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count){
	printstr();
	return PMPI_Status_set_elements_x(status, datatype, count);
}
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag, MPI_Status array_of_statuses[]){
	printstr();
	return PMPI_Testall(count, array_of_requests, flag, array_of_statuses);
}
int MPI_Testany(int count, MPI_Request array_of_requests[], int *index, int *flag, MPI_Status *status){
	printstr();
	return PMPI_Testany(count, array_of_requests, index, flag, status);
}
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status){
	printstr();
	return PMPI_Test(request, flag, status);
}
int MPI_Test_cancelled(const MPI_Status *status, int *flag){
	printstr();
	return PMPI_Test_cancelled(status, flag);
}
int MPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
	printstr();
	return PMPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
}
int MPI_Topo_test(MPI_Comm comm, int *status){
	printstr();
	return PMPI_Topo_test(comm, status);
}
int MPI_Type_commit(MPI_Datatype *type){
	printstr();
	return PMPI_Type_commit(type);
}
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_contiguous(count, oldtype, newtype);
}
int MPI_Type_create_darray(int size, int rank, int ndims, const int gsize_array[], const int distrib_array[], const int darg_array[], const int psize_array[], int order, MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_darray(size, rank, ndims, gsize_array, distrib_array, darg_array, psize_array, order, oldtype, newtype);
}
int MPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_f90_complex(p, r, newtype);
}
int MPI_Type_create_f90_integer(int r, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_f90_integer(r, newtype);
}
int MPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_f90_real(p, r, newtype);
}
int MPI_Type_create_hindexed_block(int count, int blocklength, const MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_hindexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
}
int MPI_Type_create_hindexed(int count, const int array_of_blocklengths[], const MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
}
int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
}
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn, MPI_Type_delete_attr_function *type_delete_attr_fn, int *type_keyval, void *extra_state){
	printstr();
	return PMPI_Type_create_keyval(type_copy_attr_fn, type_delete_attr_fn, type_keyval, extra_state);
}
int MPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_indexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
}
int MPI_Type_create_struct(int count, const int array_of_block_lengths[], const MPI_Aint array_of_displacements[], const MPI_Datatype array_of_types[], MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_struct(count, array_of_block_lengths, array_of_displacements, array_of_types, newtype);
}
int MPI_Type_create_subarray(int ndims, const int size_array[], const int subsize_array[], const int start_array[], int order, MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_subarray(ndims, size_array, subsize_array, start_array, order, oldtype, newtype);
}
int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_create_resized(oldtype, lb, extent, newtype);
}
int MPI_Type_delete_attr(MPI_Datatype type, int type_keyval){
	printstr();
	return PMPI_Type_delete_attr(type, type_keyval);
}
int MPI_Type_dup(MPI_Datatype type, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_dup(type, newtype);
}
int MPI_Type_free(MPI_Datatype *type){
	printstr();
	return PMPI_Type_free(type);
}
int MPI_Type_free_keyval(int *type_keyval){
	printstr();
	return PMPI_Type_free_keyval(type_keyval);
}
int MPI_Type_get_attr(MPI_Datatype type, int type_keyval, void *attribute_val, int *flag){
	printstr();
	return PMPI_Type_get_attr(type, type_keyval, attribute_val, flag);
}
int MPI_Type_get_contents(MPI_Datatype mtype, int max_integers, int max_addresses, int max_datatypes, int array_of_integers[], MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]){
	printstr();
	return PMPI_Type_get_contents(mtype, max_integers, max_addresses, max_datatypes, array_of_integers, array_of_addresses, array_of_datatypes);
}
int MPI_Type_get_envelope(MPI_Datatype type, int *num_integers, int *num_addresses, int *num_datatypes, int *combiner){
	printstr();
	return PMPI_Type_get_envelope(type, num_integers, num_addresses, num_datatypes, combiner);
}
int MPI_Type_get_extent(MPI_Datatype type, MPI_Aint *lb, MPI_Aint *extent){
	printstr();
	return PMPI_Type_get_extent(type, lb, extent);
}
int MPI_Type_get_extent_x(MPI_Datatype type, MPI_Count *lb, MPI_Count *extent){
	printstr();
	return PMPI_Type_get_extent_x(type, lb, extent);
}
int MPI_Type_get_name(MPI_Datatype type, char *type_name, int *resultlen){
	printstr();
	return PMPI_Type_get_name(type, type_name, resultlen);
}
int MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent){
	printstr();
	return PMPI_Type_get_true_extent(datatype, true_lb, true_extent);
}
int MPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *true_lb, MPI_Count *true_extent){
	printstr();
	return PMPI_Type_get_true_extent_x(datatype, true_lb, true_extent);
}
int MPI_Type_indexed(int count, const int array_of_blocklengths[], const int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_indexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
}
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype *type){
	printstr();
	return PMPI_Type_match_size(typeclass, size, type);
}
int MPI_Type_set_attr(MPI_Datatype type, int type_keyval, void *attr_val){
	printstr();
	return PMPI_Type_set_attr(type, type_keyval, attr_val);
}
int MPI_Type_set_name(MPI_Datatype type, const char *type_name){
	printstr();
	return PMPI_Type_set_name(type, type_name);
}
int MPI_Type_size(MPI_Datatype type, int *size){
	printstr();
	return PMPI_Type_size(type, size);
}
int MPI_Type_size_x(MPI_Datatype type, MPI_Count *size){
	printstr();
	return PMPI_Type_size_x(type, size);
}
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype){
	printstr();
	return PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);
}
int MPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm){
	printstr();
	return PMPI_Unpack(inbuf, insize, position, outbuf, outcount, datatype, comm);
}
int MPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name){
	printstr();
	return PMPI_Unpublish_name(service_name, info, port_name);
}
int MPI_Unpack_external (const char datarep[], const void *inbuf, MPI_Aint insize, MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype){
	printstr();
	return PMPI_Unpack_external (datarep, inbuf, insize, position, outbuf, outcount, datatype);
}
int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses){
	printstr();
	return PMPI_Waitall(count, array_of_requests, array_of_statuses);
}
int MPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status){
	printstr();
	return PMPI_Waitany(count, array_of_requests, index, status);
}
int MPI_Wait(MPI_Request *request, MPI_Status *status){
	printstr();
	return PMPI_Wait(request, status);
}
int MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
	printstr();
	return PMPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
}
