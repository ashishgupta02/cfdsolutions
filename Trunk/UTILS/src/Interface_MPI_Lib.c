/*******************************************************************************
 * File:        Interface_MPI_Lib.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>

#include "Utils.h"
#include "Check_Malloc.h"
#include "Interface_MPI_Set.h"
#include "Warn_Error.h"
#include "Interface_MPI_Lib.h"

/*******************************************************************************
 * unique data key to differentiate messages from parallel ndm
 ******************************************************************************/
static int use_mpi;

#ifdef HAVE_MPI

static int tag = 47;
static int nrequest = 0;

static MPI_Request *request = NULL;
static MPI_Status *status = NULL;
static MPI_Comm comm = MPI_COMM_WORLD;

#endif

/*******************************************************************************
 * init parallel mpi lib
 ******************************************************************************/
int lib_parallel_init(int ndomains, int mpi_mode) {
    use_mpi = mpi_mode;

    if (use_mpi) {
#ifdef HAVE_MPI
        nrequest = 0;
        request = (MPI_Request *) check_malloc(2 * ndomains * sizeof (MPI_Request));
        status = (MPI_Status *) check_malloc(2 * ndomains * sizeof (MPI_Status));
#endif
    }

    /*--------------------------------------------------------------------------
    | Note: using ndm as python extension for FlowSimulator comm is needed also
    |       if mpi_mode == FALSE (i.e.: ndomains = 1)
    --------------------------------------------------------------------------*/
#ifdef HAVE_MPI
    comm = *get_mpi_communicator();
#endif

    return use_mpi;

} /** lib_parallel_init() **/

/*******************************************************************************
 * free parallel mpi lib
 ******************************************************************************/
void lib_parallel_end(void) {

#ifdef HAVE_MPI
    nrequest = 0;
    check_free(request);
    check_free(status);
#endif

} /** lib_parallel_end() **/

/*******************************************************************************
 * MPI functions for blocking communication
 ******************************************************************************/
int mpi_send(void *sbuf, int cnt, MPI_Datatype type, int dest) {
    int err = 0;

    if (use_mpi > 0) {
#ifdef HAVE_MPI
        err = MPI_Send(sbuf, cnt, type, dest, tag, comm);
#endif
    }
    return err;
} /** mpi_send() **/

/*******************************************************************************
 * MPI functions for non blocking communication
 ******************************************************************************/
int mpi_irecv(void *rbuf, int cnt, MPI_Datatype type, int dest) {
    int err = 0;

    if (use_mpi > 0) {
#ifdef HAVE_MPI
        err = MPI_Irecv(rbuf, cnt, type, dest, tag, comm, request + nrequest++);
#endif
    }
    return err;
} /** mpi_irecv() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_isend(void *sbuf, int cnt, MPI_Datatype type, int dest) {
    int err = 0;

    if (use_mpi > 0) {
#ifdef HAVE_MPI
        err = MPI_Isend(sbuf, cnt, type, dest, tag, comm, request + nrequest++);
#endif
    }
    return err;
} /** mpi_isend() **/

/*******************************************************************************
 *
 ******************************************************************************/
void mpi_waitall(void) {
    if (use_mpi) {
#ifdef HAVE_MPI
        MPI_Waitall(nrequest, request, status);
        nrequest = 0;
#endif
    }
} /** mpi_waitall() **/

/*******************************************************************************
 * MPI functions for user defined datatypes and packing
 ******************************************************************************/
int mpi_address(void *location, MPI_Aint *address) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Address(location, address);
#endif
    }
    return err;
} /** mpi_address() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_type_commit(MPI_Datatype *datatype) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Type_commit(datatype);
#endif
    }
    return err;
} /** mpi_type_commit() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_type_free(MPI_Datatype *datatype) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Type_free(datatype);
#endif
    }
    return err;
} /** mpi_type_free() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_type_hvector(int cnt, int blocklength, MPI_Aint stride,
        MPI_Datatype oldtype, MPI_Datatype *newtype) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Type_hvector(cnt, blocklength, stride, oldtype, newtype);
#endif
    }
    return err;
} /** mpi_type_hvector() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_type_indexed(int cnt, int *blen, int *disp,
        MPI_Datatype oldtype, MPI_Datatype *newtype) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Type_indexed(cnt, blen, disp, oldtype, newtype);
#endif
    }
    return err;
} /** mpi_type_indexed() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_type_hindexed(int cnt, int *blen, MPI_Aint *disp,
        MPI_Datatype oldtype, MPI_Datatype *newtype) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Type_hindexed(cnt, blen, disp, oldtype, newtype);
#endif
    }
    return err;
} /** mpi_type_hindexed() **/

/*******************************************************************************
 *
 ******************************************************************************/
int mpi_type_struct(int cnt, int *blen, MPI_Aint *disp,
        MPI_Datatype *type, MPI_Datatype *newtype) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Type_struct(cnt, blen, disp, type, newtype);
#endif
    }
    return err;
} /** mpi_type_struct() **/

/*******************************************************************************
 * MPI functions for collective communications
 ******************************************************************************/
int mpi_allreduce(void *sbuf, void *rbuf,
        int cnt, MPI_Datatype type, MPI_Op op) {
    int err = 0;

    if (use_mpi) {
#ifdef HAVE_MPI
        err = MPI_Allreduce(sbuf, rbuf, cnt, type, op, comm);
#endif
    }
    return err;
} /** mpi_allreduce() **/

