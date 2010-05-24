/*******************************************************************************
 * File:        Interface_MPI.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "Check_Malloc.h"
#include "Interface_MPI_Lib.h"
#include "Interface_MPI_Set.h"
#include "Stopwatch.h"
#include "NDM_TypeDefs.h"
#include "Utils.h"
#include "Warn_Error.h"
#include "Interface_MPI.h"

/*******************************************************************************
 *
 *******************************************************************************/
#define DATAKEY 42
#define SBUF_SZ 128

static size_t sbuffsize = 0;
static size_t rbuffsize = 0;
static void *sbuff = NULL;
static void *rbuff = NULL;
static int use_mpi = FALSE;

static int thisdomain = 0; /*sequential default */
static int nalldomains = 1; /*sequential default */

#ifdef HAVE_MPI
static MPI_Request *request = NULL;
static MPI_Status *status = NULL;
static int reqdom = 0;
static int reqnum = 0;

static NDMDouble num_send = 0;
static NDMDouble time_send = 0;
static NDMDouble num_recv = 0;
static NDMDouble time_recv = 0;
static NDMDouble num_reduce = 0;
static NDMDouble time_reduce = 0;
#endif

static Stopwatch my_watch;

#ifdef HAVE_MPI
static MPI_Comm *my_communicator = NULL;
#endif

#ifdef HAVE_MPI
static NDMDouble dbl_recv_buffer[SBUF_SZ];
static NDMDouble dbl_send_buffer[SBUF_SZ];
static int int_send_buffer[SBUF_SZ];
static int int_recv_buffer[SBUF_SZ];
#endif

/*******************************************************************************
 *
 *******************************************************************************/
static void my_parallel_free(void);

#ifdef HAVE_MPI
static void myparallel_allreduce_int(int *v,
        int dimv,
        int comm,
        MPI_Datatype datatype,
        MPI_Op operand);
static void myparallel_allreduce_dbl(NDMDouble *v,
        int dimv,
        int comm,
        MPI_Datatype datatype,
        MPI_Op operand);

#endif

/*******************************************************************************
 * allow to set the number of processes e.g. to run a sequential prog on
 * process 0 of N processes and avoid dead locks from communication routines
 *******************************************************************************/
void set_mpi_mode(int nprocesses) {
    if (nprocesses <= 1)
        use_mpi = FALSE;
    else
        use_mpi = TRUE;

    nalldomains = nprocesses;
}

/*******************************************************************************
 *
 *******************************************************************************/
int need_mpi(void) {
    return use_mpi;
}

/*******************************************************************************
 *
 *******************************************************************************/
int mpi_thisdomain(void) {
    return thisdomain;
}

/*******************************************************************************
 * interface with naming-convention for NDM-functions in python scripts
 *******************************************************************************/
int ndm_mpi_rank(void) {
    return thisdomain;
}

/*******************************************************************************
 *
 *******************************************************************************/
int mpi_ndomains(void) {
    return nalldomains;
}

/*******************************************************************************
 * interface with naming-convention for NDM-functions in python scripts
 *******************************************************************************/
int ndm_mpi_nranks(void) {
    return nalldomains;
}

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_parallel_init(int *argc, char ***argv, int mode,
        int *ndomains, int *myid, int *mpi_mode) {
    *mpi_mode = FALSE;
    *myid = 0;
    *ndomains = 0;

    /*----------------------------------------------------------------------------
    | we assume sequential mode  <= arguments than 'mode'
    ----------------------------------------------------------------------------*/
    if (*argc <= mode) {
        ndm_msg("start sequential mode\n");

        *mpi_mode = FALSE;
        *myid = 0;
        *ndomains = 1;
    }        /*----------------------------------------------------------------------------
  | we assume mpi-mode for more than 4 commandline arguments
  ----------------------------------------------------------------------------*/
    else {
        ndm_msg("start MPI processes for parallel mode\n");

#ifdef HAVE_MPI
        my_communicator = get_mpi_communicator();

        MPI_Init(argc, argv);
        MPI_Comm_size(*my_communicator, ndomains);
        MPI_Comm_rank(*my_communicator, myid);
        *mpi_mode = TRUE;
#else
        fatal_error("MPI is not linked, missing -DHAVE_MPI, recompile!");
#endif
    }
} /** mpi_parallel_init() **/

/*******************************************************************************
 *
 *******************************************************************************/
int myparallel_init(int *argc, char ***argv, int mode) {
    int ndomains = 1;
    int this_domain = 0;
    int mpi_mode = FALSE;

    /* init MPI */
    mpi_parallel_init(argc, argv, mode, &ndomains, &this_domain, &mpi_mode);

    /* Init NDM static variables for MPI usage */
    ndm_parallel_init(ndomains, this_domain, mpi_mode);

    ndm_parallel_sync();
    return mpi_mode;
} /** myparallel_init() **/

/*******************************************************************************
 * end in case of myparallel_init had been used
 *******************************************************************************/
void myparallel_end(void) {
    ndm_parallel_end();
    mpi_parallel_end();
} /** myparallel_end() **/

/*******************************************************************************
 *
 *******************************************************************************/
void ndm_parallel_init_procs(int mpi_mode) {
    int ndomains = 1, this_domain = 0;

    if (mpi_mode) {
#ifdef HAVE_MPI
        my_communicator = get_mpi_communicator();

        MPI_Comm_size(*my_communicator, &ndomains);
        MPI_Comm_rank(*my_communicator, &this_domain);
#endif
    }
    ndm_parallel_init(ndomains, this_domain, mpi_mode);
} /** ndm_parallel_init_procs() **/

/*******************************************************************************
 *
 *******************************************************************************/
void ndm_parallel_init(int ndomains, int mydomain, int mpi_mode) {
    /*----------------------------------------------------------------------------
    | set communicator if not yet done
    ----------------------------------------------------------------------------*/
#ifdef HAVE_MPI
    my_communicator = get_mpi_communicator();
#endif

    /*----------------------------------------------------------------------------
    | reset all
    ----------------------------------------------------------------------------*/
    my_parallel_free();

    /*----------------------------------------------------------------------------
    | new init
    ----------------------------------------------------------------------------*/
    thisdomain = mydomain;
    nalldomains = ndomains;
    use_mpi = mpi_mode;

    init_stopwatch(&my_watch);
    cont_stopwatch(&my_watch);

    lib_parallel_init(ndomains, mpi_mode);

    /*----------------------------------------------------------------------------
    | check
    ----------------------------------------------------------------------------*/
    if (ndomains > 1 && !mpi_mode)
        fatal_error("Wrong initialisation: %d domains but no mpi-mode!",
            ndomains);
    if (ndomains < 0 || mydomain < 0 || mydomain >= ndomains)
        fatal_error("Wrong initialisation: domain id out of range");
} /** ndm_parallel_init() **/

/*******************************************************************************
 * if using an external communicator we have to avoid MPI_Barrier after reset
 * if e.g. the python garbage collector have already deleted
 * the initialistaion-class of the communicator: check for NULL !
 *******************************************************************************/
void ndm_parallel_sync(void) {
    if (use_mpi) {
#ifdef HAVE_MPI
        if (my_communicator != NULL && *my_communicator != MPI_COMM_NULL)
            MPI_Barrier(*my_communicator);
#endif
    }
} /** ndm_parallel_sync() **/

/*******************************************************************************
 *
 *******************************************************************************/
void ndm_parallel_print_comlog(void) {
    if (use_mpi) {
#ifdef HAVE_MPI
        if (num_send > 0 || num_recv > 0 || num_reduce > 0) {

            char p_name[MPI_MAX_PROCESSOR_NAME];
            int tmp;

            MPI_Get_processor_name(p_name, &tmp);
            MPI_Comm_rank(*my_communicator, &tmp);

            printshortline;
            ndm_msg("Communication log for domain %d on %s\n", tmp, p_name);
            ndm_msg("Wallclock runtime %.2f sec\n", read_stopwatch(&my_watch));
            ndm_msg("%g sends in %.2g sec,  ",
                    num_send, time_send);
            ndm_msg("%g receives in %.2g sec,  ",
                    num_recv, time_recv);
            ndm_msg("%g reduces in %.2g sec\n",
                    num_reduce, time_reduce);

            printshortline;
        }
#endif

    }
} /** ndm_parallel_print_comlog() **/

/*******************************************************************************
 *
 *******************************************************************************/
void ndm_parallel_end(void) {
    if (!use_mpi)
        ndm_msg("Wallclock runtime %.2f sec\n", read_stopwatch(&my_watch));
    else {
        ndm_parallel_sync();
        ndm_msg("MPI_Finalize\n");
        ndm_parallel_print_comlog();
        my_parallel_free();
        lib_parallel_end();
    }

    stop_stopwatch(&my_watch);
} /** ndm_parallel_end() **/

/*******************************************************************************
 *
 *******************************************************************************/
void ndm_parallel_abort(void) {
    if (use_mpi) {
#ifdef HAVE_MPI
        const int myid = mpi_thisdomain();
        const int ndom = mpi_ndomains();
        int hard_abort = TRUE, flag, i;
        MPI_Request req;
        MPI_Status st;

        if (myid == 0) {
            MPI_Isend(&i, 1, MPI_INT, (myid + 1) % ndom, 2706, MPI_COMM_WORLD, &req);
            for (flag = 0, i = 0; !flag && i < 10000; i++) {
                msleep(1);
                MPI_Iprobe((myid - 1 + ndom) % ndom, 2706, MPI_COMM_WORLD, &flag, &st);
            }
            if (flag) {
                MPI_Recv(&i, 1, MPI_INT, (myid - 1 + ndom) % ndom, 2706, MPI_COMM_WORLD, &st);
                hard_abort = FALSE;
            }
        } else { /* i.e., myid > 0  */
            for (flag = 0, i = 0; !flag && i < 10000; i++) {
                msleep(1);
                MPI_Iprobe((myid - 1 + ndom) % ndom, 2706, MPI_COMM_WORLD, &flag, &st);
            }
            if (flag) {
                MPI_Recv(&i, 1, MPI_INT, (myid - 1 + ndom) % ndom, 2706, MPI_COMM_WORLD, &st);
                MPI_Send(&i, 1, MPI_INT, (myid + 1) % ndom, 2706, MPI_COMM_WORLD);
                hard_abort = FALSE;
            }
        }
        if (hard_abort) {
            ndm_print(stderr, 0, "\nABORT on domain (MPI rank) %d\n\n", myid);
            fflush(stderr);
            ndm_print_close_stdout_stderr();
            MPI_Abort(MPI_COMM_WORLD, 4711);
        } else {
            if (myid == 0) ndm_print(stderr, 0, "\nABORT on all domains\n\n");
            fflush(stderr);
            ndm_parallel_end();
            mpi_parallel_end();
            ndm_print_close_stdout_stderr();
        }
#endif
    }
    exit(EXIT_FAILURE);
} /** ndm_parallel_abort() **/

/*******************************************************************************
 * reset signal handler for mpi_parallel_end(void)
 * ... we found a case where sig 10 was send aound when finalizing mpi
 *     obviously due to previous errors. In order to avoid misleading
 *     signal-messages we reset the signal handler before terminating
 *******************************************************************************/
static void term_signal_handler(int sig) {

#ifdef HAVE_MPI
    int tmp;
    MPI_Comm_rank(*my_communicator, &tmp);
    ndm_errmsg("Signal %d before/from MPI_Finalize() on process %d\n",
            sig, thisdomain);
    ndm_msg("Signal %d before/from MPI_Finalize() on process %d\n",
            sig, tmp);
#else
    return;
#endif
} /** term_signal_handler() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_parallel_end(void) {
    if (use_mpi) {

        /*--------------------------------------------------------------------------
        | Barrier needed for IBM in Corbeil with  POE MPI
        | Signal handling and MPI_Finalize() can otherwise lead to hanging processes
        --------------------------------------------------------------------------*/
#ifdef HAVE_MPI
        if (my_communicator != NULL && *my_communicator != MPI_COMM_NULL) {
            MPI_Barrier(*my_communicator);
        }
#endif

        /*--------------------------------------------------------------------------
        | reset signal handler messages
        --------------------------------------------------------------------------*/
        signal(SIGUSR1, term_signal_handler);
        signal(SIGUSR2, term_signal_handler);

        /*--------------------------------------------------------------------------
        | finalize mpi
        --------------------------------------------------------------------------*/
#ifdef HAVE_MPI
        if (my_communicator != NULL && *my_communicator != MPI_COMM_NULL) {
            MPI_Barrier(*my_communicator);
            MPI_Finalize();
            my_communicator = NULL;
        }
#endif
    }
    use_mpi = FALSE;
} /** mpi_parallel_end() **/

/*******************************************************************************
 * free/reset ststic data
 *******************************************************************************/
static void my_parallel_free(void) {
    rbuffsize = 0;
    sbuffsize = 0;
    check_free(sbuff);
    check_free(rbuff);

    thisdomain = 0;
    nalldomains = 1;

#ifdef HAVE_MPI
    check_free(status);
    check_free(request);
    reqnum = 0;
    reqdom = 0;
    num_send = 0;
    num_recv = 0;
    time_send = 0;
    time_recv = 0;
#endif
} /** my_parallel_free() **/

/*******************************************************************************
 *
 *******************************************************************************/
void *my_send_commbuffer(size_t size) {
    void *buff = sbuff;

    if (use_mpi) {
        if (size > sbuffsize) {
            check_free(sbuff);
            sbuffsize = size;
            sbuff = check_malloc(sbuffsize);
            buff = sbuff;
        }
    } else {
        if (size > 0) {
            terminate("No communication buffers for dummy communication!");
        }
    }

    return buff;
} /** my_send_commbuffer() **/

/*******************************************************************************
 *
 *******************************************************************************/
void *my_recv_commbuffer(size_t size) {
    void *buff = rbuff;

    if (use_mpi) {
        if (size > rbuffsize) {
            check_free(rbuff);
            rbuffsize = size;
            rbuff = check_malloc(rbuffsize);
            buff = rbuff;
        }
    } else {
        if (size > 0) {
            terminate("No communication buffers for dummy communication!");
        }
    }

    return buff;
} /** my_recv_commbuffer() **/

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_send(void *buf, size_t size, int to, int key) {
    if (size == 0)
        return;

    if (use_mpi) {
#ifdef HAVE_MPI
        NDMDouble starttime = MPI_Wtime();

        MPI_Send(buf, (int) size, MPI_BYTE, to, key, *my_communicator);

        num_send++;
        time_send += MPI_Wtime() - starttime;
#endif
    } else {
        terminate("No send for dummy communication!");
    }
} /** myparallel_send() **/

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_recv(void *buf, size_t size, int from, int key) {

    if (size == 0)
        return;

    if (use_mpi) {
#ifdef HAVE_MPI
        NDMDouble starttime = MPI_Wtime();
        MPI_Status stat;

        MPI_Recv(buf, (int) size, MPI_BYTE, from, key, *my_communicator, &stat);

        num_recv++;
        time_recv += MPI_Wtime() - starttime;
#endif
    } else {
        terminate("No receive for dummy communication!");
    }
} /** myparallel_recv() **/

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_int_globalsum(int *v, int dimv, int key) {
#ifdef HAVE_MPI
    myparallel_allreduce_int(v, dimv, -1, MPI_INT, MPI_SUM);
#endif
}

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_dbl_globalsum(NDMDouble *v, int dimv, int key) {
#ifdef HAVE_MPI
    myparallel_allreduce_dbl(v, dimv, -1, MPI_DOUBLE, MPI_SUM);
#endif
}

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_int_globalmax(int *v, int dimv, int key) {
#ifdef HAVE_MPI
    myparallel_allreduce_int(v, dimv, -1, MPI_INT, MPI_MAX);
#endif
}

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_dbl_globalmax(NDMDouble *v, int dimv, int key) {
#ifdef HAVE_MPI
    myparallel_allreduce_dbl(v, dimv, -1, MPI_DOUBLE, MPI_MAX);
#endif
}

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_int_globalmin(int *v, int dimv, int key) {
#ifdef HAVE_MPI
    myparallel_allreduce_int(v, dimv, -1, MPI_INT, MPI_MIN);
#endif
}

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_dbl_globalmin(NDMDouble *v, int dimv, int key) {
#ifdef HAVE_MPI
    myparallel_allreduce_dbl(v, dimv, -1, MPI_DOUBLE, MPI_MIN);
#endif
}

/*******************************************************************************
 * collect sum of integer variables on domains of a splitted communicator
 *******************************************************************************/
void myparallel_int_localsum(int *v, int dimv, int comm) {
#ifdef HAVE_MPI
    myparallel_allreduce_int(v, dimv, comm, MPI_INT, MPI_SUM);
#endif
}

/*******************************************************************************
 * collect sum of double variables on domains of a splitted communicator
 *******************************************************************************/
void myparallel_dbl_localsum(NDMDouble *v, int dimv, int comm) {
#ifdef HAVE_MPI
    myparallel_allreduce_dbl(v, dimv, comm, MPI_DOUBLE, MPI_SUM);
#endif
}

/*******************************************************************************
 * compute maxima of integer variables on domains of a splitted communicator
 *******************************************************************************/
void myparallel_int_localmax(int *v, int dimv, int comm) {
#ifdef HAVE_MPI
    myparallel_allreduce_int(v, dimv, comm, MPI_INT, MPI_MAX);
#endif
}

/*******************************************************************************
 * compute maxima of double variables on domains of a splitted communicator
 *******************************************************************************/
void myparallel_dbl_localmax(NDMDouble *v, int dimv, int comm) {
#ifdef HAVE_MPI
    myparallel_allreduce_dbl(v, dimv, comm, MPI_DOUBLE, MPI_MAX);
#endif
}

/*******************************************************************************
 * compute minima of integer variables on domains of a splitted communicator
 *******************************************************************************/
void myparallel_int_localmin(int *v, int dimv, int comm) {
#ifdef HAVE_MPI
    myparallel_allreduce_int(v, dimv, comm, MPI_INT, MPI_MIN);
#endif
}

/*******************************************************************************
 * compute minima of double variables on domains of a splitted communicator
 *******************************************************************************/
void myparallel_dbl_localmin(NDMDouble *v, int dimv, int comm) {
#ifdef HAVE_MPI
    myparallel_allreduce_dbl(v, dimv, comm, MPI_DOUBLE, MPI_MIN);
#endif
}

/*******************************************************************************
 * Boolean OR over domains.
 *******************************************************************************/
void exchange_boolean_or(int *flag) {
    int vglobal[1];
    vglobal[0] = (*flag) ? 1 : 0;
    myparallel_int_globalmax(vglobal, 1, 4321);
    *flag = (vglobal[0] > 0) ? TRUE : FALSE;
}

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_send_string(char *buf, int dest, int tag) {
    if (use_mpi) {
#ifdef HAVE_MPI

        int len = buf != NULL ? strlen(buf) + 1 : 0;

        if (len > 0)
            MPI_Send(buf, len, MPI_CHAR, dest, tag, *my_communicator);
        else
            fatal_error("empty string");

#endif
    }
} /** myparallel_send_string() **/

/*******************************************************************************
 *
 *******************************************************************************/
void myparallel_recv_string(char *buf, int size, int tag) {
    if (use_mpi) {
#ifdef HAVE_MPI
        MPI_Status my_status;

        MPI_Recv(buf, size, MPI_CHAR, thisdomain, tag,
                *my_communicator, &my_status);
#endif
    }
} /** myparallel_recv_string() **/

/*******************************************************************************
 *
 *******************************************************************************/
void my_nb_sendrecv_init(int ndomains) {
    if (use_mpi) {
#ifdef HAVE_MPI

        reqnum = 0;

        if (ndomains > 0 && ndomains != reqdom) {
            size_t sz_r = 2 * ndomains * sizeof (MPI_Request);
            size_t sz_s = 2 * ndomains * sizeof (MPI_Status);

            check_free(request);
            request = (MPI_Request *) check_malloc(sz_r);

            check_free(status);
            status = (MPI_Status *) check_malloc(sz_s);
            reqdom = ndomains;
        }
#endif
    }
} /** my_nb_sendrecv_init() **/

/*******************************************************************************
 *
 *******************************************************************************/
void my_nb_recv(void *rbuf, size_t size, int dom) {
    if (use_mpi && size > 0) {
#ifdef HAVE_MPI
        NDMDouble starttime = MPI_Wtime();

        MPI_Irecv(rbuf, (int) size, MPI_BYTE, dom, DATAKEY, *my_communicator,
                (&request[reqnum++]));

        num_recv++;
        time_recv += MPI_Wtime() - starttime;
#endif
    }
} /** my_nb_recv() **/

/*******************************************************************************
 *
 *******************************************************************************/
void my_nb_send(void *sbuf, size_t size, int dom) {
    if (use_mpi && size > 0) {
#ifdef HAVE_MPI
        NDMDouble starttime = MPI_Wtime();

        MPI_Isend(sbuf, (int) size, MPI_BYTE, dom, DATAKEY, *my_communicator,
                (&request[reqnum++]));

        num_send++;
        time_send += MPI_Wtime() - starttime;
#endif
    }
} /** my_nb_send() **/

/*******************************************************************************
 *
 *******************************************************************************/
void my_nb_waitall(void) {
    if (use_mpi) {
#ifdef HAVE_MPI
        MPI_Waitall(reqnum, request, status);
#endif
    }
} /** my_nb_waitall() **/

/*******************************************************************************
 * compute the sum/max/min of intergers of all domains of a given communicator
 *******************************************************************************/
#ifdef HAVE_MPI

static void myparallel_allreduce_int(int *v,
        int dimv,
        int comm,
        MPI_Datatype datatype,
        MPI_Op operand) {
    if (use_mpi) {
        MPI_Comm *this_communicator = NULL;
        NDMDouble starttime;

        int i;
        int *obuf = NULL;
        int *sbuf = NULL;

        /*--------------------------------------------------------------------------
        | set the communicator
        --------------------------------------------------------------------------*/
        if (comm == -1)
            this_communicator = get_mpi_communicator();
        else
            this_communicator = get_mpi_split_communicator(comm);

        if (*this_communicator == MPI_COMM_NULL)
            terminate("Invalid communicator for this domain");

        /*--------------------------------------------------------------------------
        | copy data to buffer, create a bigger buffer if dimv exceeds SBUF_SZ
        --------------------------------------------------------------------------*/
        if (dimv > SBUF_SZ) {
            obuf = ndmint1_mem(dimv);
            sbuf = ndmint1_mem(dimv);
        } else {
            obuf = int_recv_buffer;
            sbuf = int_send_buffer;
        }

        for (i = 0; i < dimv; i++)
            sbuf[i] = v[i];

        /*--------------------------------------------------------------------------
        | count time and call Allreduce
        --------------------------------------------------------------------------*/
        starttime = MPI_Wtime();

        MPI_Allreduce(sbuf, obuf, dimv, datatype, operand, *this_communicator);

        num_reduce++;
        time_reduce += MPI_Wtime() - starttime;

        /*--------------------------------------------------------------------------
        | copy data back
        --------------------------------------------------------------------------*/
        for (i = 0; i < dimv; i++)
            v[i] = obuf[i];

        if (dimv > SBUF_SZ) {
            check_free(obuf);
            check_free(sbuf);
        }
    } /** if(use_mpi) **/

} /** myparallel_allreduce_int() **/
#endif

/*******************************************************************************
 * compute the sum/max/min of doubles of all domains of a given communicator
 *******************************************************************************/
#ifdef HAVE_MPI

static void myparallel_allreduce_dbl(NDMDouble *v,
        int dimv,
        int comm,
        MPI_Datatype datatype,
        MPI_Op operand) {
    if (use_mpi) {
        MPI_Comm *this_communicator = NULL;
        NDMDouble starttime;

        int i;
        NDMDouble *obuf = NULL;
        NDMDouble *sbuf = NULL;

        /*--------------------------------------------------------------------------
        | set the communicator
        --------------------------------------------------------------------------*/
        if (comm == -1)
            this_communicator = get_mpi_communicator();
        else
            this_communicator = get_mpi_split_communicator(comm);

        if (*this_communicator == MPI_COMM_NULL)
            terminate("Invalid communicator for this domain");

        /*--------------------------------------------------------------------------
        | copy data to buffer, create a bigger buffer if dimv exceeds SBUF_SZ
        --------------------------------------------------------------------------*/
        if (dimv > SBUF_SZ) {
            obuf = ndmdouble1_mem(dimv);
            sbuf = ndmdouble1_mem(dimv);
        } else {
            obuf = dbl_recv_buffer;
            sbuf = dbl_send_buffer;
        }

        for (i = 0; i < dimv; i++)
            sbuf[i] = v[i];

        /*--------------------------------------------------------------------------
        | count time and call Allreduce
        --------------------------------------------------------------------------*/
        starttime = MPI_Wtime();

        MPI_Allreduce(sbuf, obuf, dimv, datatype, operand, *this_communicator);

        num_reduce++;
        time_reduce += MPI_Wtime() - starttime;

        /*--------------------------------------------------------------------------
        | copy data back
        --------------------------------------------------------------------------*/
        for (i = 0; i < dimv; i++)
            v[i] = obuf[i];

        if (dimv > SBUF_SZ) {
            check_free(obuf);
            check_free(sbuf);
        }
    } /** if(use_mpi) **/

} /** myparallel_allreduce_dbl() **/
#endif

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_gather_data(void *send, void *recv, int dim, int data_type) {
    if (use_mpi) {
#ifdef HAVE_MPI

        if (data_type == 0)
            MPI_Gather(send, dim, MPI_INT, recv, dim, MPI_INT, 0, *my_communicator);

        if (data_type == 1)
            MPI_Gather(send, dim, MPI_DOUBLE, recv, dim, MPI_DOUBLE, 0,
                *my_communicator);

#endif
    }

} /** mpi_gather_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_gatherv_data(void *send, void *recv, int sdim, int *rdim, int *disp,
        int data_type) {
    if (use_mpi) {
#ifdef HAVE_MPI
        if (data_type == 0)
            MPI_Gatherv(send, sdim, MPI_INT, recv, rdim, disp, MPI_INT, 0,
                *my_communicator);

        if (data_type == 1)
            MPI_Gatherv(send, sdim, MPI_DOUBLE, recv, rdim, disp, MPI_DOUBLE, 0,
                *my_communicator);
#endif
    }

} /** mpi_gatherv_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_scatter_data(void *send, void *recv, int dim, int data_type) {
    if (use_mpi) {
#ifdef HAVE_MPI

        if (data_type == 0)
            MPI_Scatter(send, dim, MPI_INT, recv, dim, MPI_INT, 0,
                *my_communicator);

        if (data_type == 1)
            MPI_Scatter(send, dim, MPI_DOUBLE, recv, dim, MPI_DOUBLE, 0,
                *my_communicator);

#endif
    }

} /** mpi_scatter_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_scatterv_data(void *send, void *recv, int *sdim, int rdim, int *disp,
        int data_type) {
    if (use_mpi) {
#ifdef HAVE_MPI

        if (data_type == 0)
            MPI_Scatterv(send, sdim, disp, MPI_INT, recv, rdim, MPI_INT, 0,
                *my_communicator);

        if (data_type == 1)
            MPI_Scatterv(send, sdim, disp, MPI_DOUBLE, recv, rdim, MPI_DOUBLE, 0,
                *my_communicator);

#endif
    }

} /** mpi_scatterv_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_send_data(void *send, int dim, int to, int data_type) {

    if (use_mpi) {
#ifdef HAVE_MPI

        NDMDouble starttime = MPI_Wtime();

        if (data_type == 0)
            MPI_Send(send, dim, MPI_INT, to, DATAKEY, *my_communicator);

        if (data_type == 1)
            MPI_Send(send, dim, MPI_DOUBLE, to, DATAKEY, *my_communicator);

        num_send++;
        time_send += MPI_Wtime() - starttime;

#endif
    }

} /** mpi_send_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_recv_data(void *recv, int dim, int from, int data_type) {

    if (use_mpi) {
#ifdef HAVE_MPI

        NDMDouble starttime = MPI_Wtime();
        MPI_Status mystat;

        if (data_type == 0)
            MPI_Recv(recv, dim, MPI_INT, from, DATAKEY, *my_communicator, &mystat);

        if (data_type == 1)
            MPI_Recv(recv, dim, MPI_DOUBLE, from, DATAKEY, *my_communicator, &mystat);

        num_recv++;
        time_recv += MPI_Wtime() - starttime;

#endif
    }
} /** mpi_recv_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_allgather_data(void *send, void *recv, int dim, int data_type) {
    if (use_mpi) {
#ifdef HAVE_MPI

        if (data_type == 0) {
            MPI_Allgather(send, dim, MPI_INT, recv, dim, MPI_INT,
                    *my_communicator);
        }

        if (data_type == 1) {
            MPI_Allgather(send, dim, MPI_DOUBLE, recv, dim, MPI_DOUBLE,
                    *my_communicator);
        }

#endif
    }
} /** mpi_allgather_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_allgatherv_data(void *send, int sdim, void *recv, int* rdim,
        int data_type) {
    int ndomains = mpi_ndomains();
    int* offset = check_malloc(ndomains * sizeof (int));
    int i;

    offset[0] = 0;
    for (i = 1; i < ndomains; i++) {
        offset[i] = offset[i - 1] + rdim[i - 1];
    }

    if (use_mpi) {
#ifdef HAVE_MPI

        if (data_type == 0) {
            MPI_Allgatherv(send, sdim, MPI_INT, recv, rdim, offset,
                    MPI_INT, *my_communicator);
        }

        if (data_type == 1) {
            MPI_Allgatherv(send, sdim, MPI_DOUBLE, recv, rdim, offset,
                    MPI_DOUBLE, *my_communicator);
        }

#endif
    }
    check_free(offset);
} /** mpi_allgatherv_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_broadcast_data(void* buffer, int dim, int root, int data_type) {
    if (use_mpi) {
#ifdef HAVE_MPI

        if (data_type == 0) {
            MPI_Bcast(buffer, dim, MPI_INT, root, *my_communicator);
        }

        if (data_type == 1) {
            MPI_Bcast(buffer, dim, MPI_DOUBLE, root, *my_communicator);
        }

#endif
    }
} /** mpi_broadcast_data() **/

/*******************************************************************************
 *
 *******************************************************************************/
void mpi_localbroadcast_data(void* buffer, int dim, int root, int data_type,
        int comm) {
    if (use_mpi) {
#ifdef HAVE_MPI
        MPI_Comm *this_communicator = NULL;

        /*--------------------------------------------------------------------------
        | set the communicator
        --------------------------------------------------------------------------*/
        if (comm == -1)
            this_communicator = get_mpi_communicator();
        else
            this_communicator = get_mpi_split_communicator(comm);

        if (data_type == 0) {
            MPI_Bcast(buffer, dim, MPI_INT, root, *this_communicator);
        }

        if (data_type == 1) {
            MPI_Bcast(buffer, dim, MPI_DOUBLE, root, *this_communicator);
        }

#endif
    }
} /** mpi_localbroadcast_data() **/

