/*******************************************************************************
 * File:        Interface_MPI_Lib.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _INTERFACE_MPI_LIB_H
#define _INTERFACE_MPI_LIB_H

#include "NDM_TypeDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#else

typedef enum {
  MPI_INT,
  MPI_DOUBLE,
  MPI_UB,
  MPI_2INT
} MPI_Datatype;

typedef enum {
  MPI_BOTTOM
} MPI_Aint;

typedef enum {
  MPI_MAX,
  MPI_MIN,
  MPI_SUM,
  MPI_MAXLOC
} MPI_Op;

#endif

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif


/*******************************************************************************
* init library
*******************************************************************************/
int lib_parallel_init(int ndomains, int mpi_mode);

/*******************************************************************************
* free library
*******************************************************************************/
void lib_parallel_end(void);

/*******************************************************************************
* MPI functions for blocking communication
*******************************************************************************/
int mpi_send(void *sbuf, int cnt, MPI_Datatype type, int dest);

/*******************************************************************************
* MPI functions for non blocking communication
*******************************************************************************/
int mpi_irecv(void *rbuf, int cnt, MPI_Datatype type, int dest);

int mpi_isend(void *sbuf, int cnt, MPI_Datatype type, int dest);

void mpi_waitall(void);

/*******************************************************************************
* MPI functions for user defined datatypes and packing
*******************************************************************************/
int mpi_address(void *location, MPI_Aint *address);


int mpi_type_commit(MPI_Datatype *datatype);


int mpi_type_free(MPI_Datatype *datatype);


int mpi_type_hvector(int cnt, int blocklength, MPI_Aint stride,
                     MPI_Datatype oldtype, MPI_Datatype *newtype);

int mpi_type_indexed(int cnt, int *blen, int *disp,
                     MPI_Datatype oldtype, MPI_Datatype *newtype);

int mpi_type_hindexed(int cnt, int *blen, MPI_Aint *disp,
                      MPI_Datatype oldtype, MPI_Datatype *newtype);

int mpi_type_struct(int cnt, int *blen, MPI_Aint *disp,
                    MPI_Datatype *type, MPI_Datatype *newtype);

/*******************************************************************************
* MPI functions for collective communications
*******************************************************************************/
int mpi_allreduce(void *sbuf, void *rbuf,
                  int cnt, MPI_Datatype type, MPI_Op op);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _INTERFACE_MPI_LIB_H */
