/*******************************************************************************
 * File:        Interface_MPI_Set.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdio.h>

#include "Utils.h"
#include "Warn_Error.h"
#include "Interface_MPI_Set.h"

/*------------------------------------------------------------------------------
| the following block is active if an mpi library is linked
------------------------------------------------------------------------------*/
#ifdef HAVE_MPI

#include "mpi.h"

/*------------------------------------------------------------------------------
| static MPI variables
------------------------------------------------------------------------------*/
static MPI_Comm *my_communicator = NULL;
static MPI_Comm default_communicator = MPI_COMM_WORLD;

static MPI_Comm my_split_communicator[MAX_NO_COMM] = {MPI_COMM_NULL};
static int my_split_comm_size[MAX_NO_COMM] = {1};

/*******************************************************************************
 * set communicator at start otherwise MPI_COMM_WORLD is used
 *******************************************************************************/
void set_mpi_communicator(MPI_Comm *set_communicator) {
    my_communicator = set_communicator;
}

/*******************************************************************************
 * provide a pointer to the actual mpi communicator
 *******************************************************************************/
MPI_Comm *get_mpi_communicator(void) {
    if (my_communicator == NULL) {
        my_communicator = &default_communicator;
    }
    return (my_communicator);
}

/*******************************************************************************
 * provide a pointer to the splitted mpi communicator with the index i_comm
 *******************************************************************************/
MPI_Comm *get_mpi_split_communicator(int i_comm) {
    if (i_comm >= MAX_NO_COMM)
        terminate("The number of required splitted communicators exceeds\n"
            "the maximum defined by MAX_NO_COMM");

    return &my_split_communicator[i_comm];
}

/*------------------------------------------------------------------------------
| the following block is active if NO mpi library is linked
------------------------------------------------------------------------------*/
#else  /** HAVE_MPI **/

static int my_null = MPI_COMM_NULL;

/*******************************************************************************
 * dummy get function
 *******************************************************************************/
int *get_mpi_communicator(void) {
    return &my_null;
}

/*******************************************************************************
 * dummy get function
 *******************************************************************************/
int *get_mpi_split_communicator(int i_comm) {
    return &my_null;
}

#endif /** HAVE_MPI **/

/*******************************************************************************
 * split mpi communicators (empty function if mpi is not available)
 *******************************************************************************/
void split_communicator(int i_comm, int color, int key) {

#ifdef HAVE_MPI
    int initilized = 0;

    MPI_Initialized(&initilized);

    if (!initilized)
        return;

    if (color < 0)
        color = MPI_UNDEFINED;

    MPI_Comm_split(*my_communicator, color, key, &my_split_communicator[i_comm]);

    if (my_split_communicator[i_comm] != MPI_COMM_NULL)
        MPI_Comm_size(my_split_communicator[i_comm], &my_split_comm_size[i_comm]);

#endif

} /** split_communicator() **/

/*******************************************************************************
 * split mpi communicators (empty function if mpi is not available)
 *******************************************************************************/
int get_split_commsize(int i_comm) {
#ifdef HAVE_MPI
    return my_split_comm_size[i_comm];
#else
    return 1;
#endif
}

