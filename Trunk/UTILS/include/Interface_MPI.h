/*******************************************************************************
 * File:        Interface_MPI.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _INTERFACE_MPI_H
#define _INTERFACE_MPI_H

#include <stddef.h>

#include "NDM_TypeDefs.h"

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*******************************************************************************
*
*******************************************************************************/
#define MPIKEY_DEFAULT 4711

/*******************************************************************************
* allow to set the number of processes e.g. to run a sequential prog on
* process 0 of N processes
*******************************************************************************/
void set_mpi_mode(int nprocesses);

/*******************************************************************************
* return status if mpi is needed (MPI_init is done!)
*******************************************************************************/
int need_mpi(void);
int mpi_thisdomain(void);
int mpi_ndomains(void);

int ndm_mpi_rank(void);
int ndm_mpi_nranks(void);
/*******************************************************************************
* set number of maximum commandline arguments
*******************************************************************************/
#define STAND_ALONE_MODE   3
#define PYTHON_SCRIPT_MODE 4

/*******************************************************************************
* init/ ndm_parallel and mpi
* !!! use argument mode = STAND_ALONE_MODE or mode = PYHTHON_SCRIPT_MODE
* !!! for switch between sequential or parallel-mpi start up
*******************************************************************************/
int myparallel_init(int *argc, char ***argv, int mode);

/*******************************************************************************
* same then myparallel_init() without ndm_parallel_init()
*******************************************************************************/
void mpi_parallel_init(int *argc, char ***argv, int mode,
                       int *ndomains, int *myid, int *mpi_mode);

/*******************************************************************************
* end ndm_parallel and mpi
*******************************************************************************/
void myparallel_end(void);
void mpi_parallel_end(void); /* end mpi only */

/*******************************************************************************
* init/end ndm_parallel if mpi is initiated elsewhere
*******************************************************************************/
void ndm_parallel_init_procs(int mpi_mode);  /*ndm_parallel_init() by MPI_RANK*/
void ndm_parallel_init(int ndomains, int mydomain, int mpi_mode);
void ndm_parallel_end(void);
void ndm_parallel_abort(void);

void ndm_parallel_print_comlog(void);

/*******************************************************************************
* sync all processes (set barrier)
*******************************************************************************/
void ndm_parallel_sync(void);


/*******************************************************************************
* Get a buffer. One for send one for recv. Do not free it!
*******************************************************************************/
void *my_send_commbuffer(size_t size);
void *my_recv_commbuffer(size_t size);

/*******************************************************************************
* blocking send
* The buffer 'buf' must be obtained from 'my_send/recv_commbuffer()'
*******************************************************************************/
void myparallel_send(void *buf, size_t size, int to, int key);

/*******************************************************************************
* blocking receive
* The buffer 'buf' must be obtained from 'myparallel_communication_buffer()'
*******************************************************************************/
void myparallel_recv(void *buf, size_t size, int from, int key);

/*******************************************************************************
* non-blocking send/recv
*******************************************************************************/
void my_nb_sendrecv_init(int ndomains);
void my_nb_recv(void *rbuf, size_t size, int dom);
void my_nb_send(void *sbuf, size_t size, int dom);
void my_nb_waitall(void);

/*******************************************************************************
* functions for computing the sum/min/max of integer or double variables
* on all domains of the global communicator
*******************************************************************************/
void myparallel_int_globalsum(int       *v, int dimv, int key);
void myparallel_dbl_globalsum(NDMDouble *v, int dimv, int key);

#define myparallel_globalsum(v, dimv, key) \
        (myparallel_dbl_globalsum(v, dimv, key))

void myparallel_int_globalmax(int       *v, int dimv, int key);
void myparallel_dbl_globalmax(NDMDouble *v, int dimv, int key);

#define myparallel_globalmax(v, dimv, key) \
        (myparallel_dbl_globalmax(v, dimv, key))

void myparallel_int_globalmin(int       *v, int dimv, int key);
void myparallel_dbl_globalmin(NDMDouble *v, int dimv, int key);

#define myparallel_globalmin(v, dimv, key) \
        (myparallel_dbl_globalmin(v, dimv, key))

/*******************************************************************************
* functions for computing the sum/min/max of integer or double variables
* on all domains of a splitted (reduced) communicator
*******************************************************************************/
void myparallel_int_localsum(int       *v, int dimv, int comm);
void myparallel_dbl_localsum(NDMDouble *v, int dimv, int comm);

void myparallel_int_localmax(int       *v, int dimv, int comm);
void myparallel_dbl_localmax(NDMDouble *v, int dimv, int comm);

void myparallel_int_localmin(int       *v, int dimv, int comm);
void myparallel_dbl_localmin(NDMDouble *v, int dimv, int comm);

/*******************************************************************************
*  mpi communication functions
*******************************************************************************/
void myparallel_recv_string(char *buf, int size, int tag);
void myparallel_send_string(char *buf, int dest, int tag);

void mpi_gather_data(void *send, void *recv, int dim, int data_type);
void mpi_gatherv_data(void *send, void *recv, int sdim, int *rdim, int *disp,
                      int data_type);
void mpi_scatter_data(void *send, void *recv, int dim, int data_type);
void mpi_scatterv_data(void *send, void *recv, int *sdim, int rdim, int *disp,
                       int data_type);
void mpi_send_data(void *send, int dim, int to, int data_type);
void mpi_recv_data(void *recv, int dim, int from, int data_type);


void mpi_allgather_data(void *send, void *recv, int dim, int data_type);
void mpi_allgatherv_data(void *send, int sdim, void *recv, int* rdim,
                         int data_type);
void mpi_broadcast_data(void* buffer, int dim, int root, int data_type);
void mpi_localbroadcast_data(void* buffer, int dim, int root, int data_type, 
                             int comm);

/*******************************************************************************
* Boolean OR over domains.
*******************************************************************************/
void exchange_boolean_or(int *flag);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _INTERFACE_MPI_H */

