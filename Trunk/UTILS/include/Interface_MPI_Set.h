/*******************************************************************************
 * File:        Interface_MPI_Set.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _INTERFACE_MPI_SET_H
#define _INTERFACE_MPI_SET_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| the following block is active if an mpi library is linked
------------------------------------------------------------------------------*/
#ifdef HAVE_MPI

# include "mpi.h"

#define MAX_NO_COMM 128

/*******************************************************************************
* set communicator at start otherwise MPI_COMM_WORLD is used
*******************************************************************************/
void set_mpi_communicator(MPI_Comm *set_communicator);

/*******************************************************************************
* provide a pointer to the actual mpi communicator
*******************************************************************************/
MPI_Comm *get_mpi_communicator(void);

/*******************************************************************************
* provide a pointer to the splitted mpi communicator with the index i_comm
*******************************************************************************/
MPI_Comm *get_mpi_split_communicator(int i_comm);

/*------------------------------------------------------------------------------
| the following block is active if NO mpi library is linked
------------------------------------------------------------------------------*/
#else  /** HAVE_MPI **/

#ifndef MPI_COMM_NULL
# define MPI_COMM_NULL -1
#endif

int *get_mpi_communicator(void);

int *get_mpi_split_communicator(int i_comm);

#endif /** HAVE_MPI **/

/*******************************************************************************
* split mpi communicators (empty function if mpi is not available)
*******************************************************************************/
void split_communicator(int i_comm, int color, int key);

/*******************************************************************************
* split mpi communicators (empty function if mpi is not available)
*******************************************************************************/
int get_split_commsize(int i_comm);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /** _INTERFACE_MPI_SET_H **/

