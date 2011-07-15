/*******************************************************************************
 * File:        Warn_Error.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Check_Malloc.h"
#include "Interface_MPI.h"
#include "Logging.h"
#include "Utils.h"
#include "Log.h"
#include "Warn_Error.h"
#include "Output_Level.h"

static int error_flag = 0;

/*******************************************************************************
 *
 *******************************************************************************/
_message_fun _print_position_and_return(const char *file, int line,
        _message_fun call) {
    if (error_flag == 0 && (have_logfile() || check_errstream() != NULL)) {
        FILE *f = get_errstream();
        ndm_print(f, PRINT_ALWAYS, "\nFILE '%s', LINE %d:\n", file, line);
    }

    if (error_flag == 0) {
        FILE *f = get_outstream();
        ndm_print(f, PRINT_ALWAYS, "\nFILE '%s', LINE %d:\n", file, line);
    }
    return call;
} /** _print_position_and_return() **/

/*******************************************************************************
 * DO NOT CALL THIS TWICE FROM INBETWEEN va_start() and va_end()
 * this is not save for old ANSI-Standard!
 *******************************************************************************/
static void write_out(const char *prefix, const char *fmt, va_list arg) {
    FILE *f = get_outstream();

    ndm_print(f, PRINT_ALWAYS, prefix);
    ndm_printlist(f, PRINT_ALWAYS, fmt, arg);
    ndm_print(f, PRINT_ALWAYS, "\n\n");
    fflush(f);
} /** write_out() **/

/*******************************************************************************
 * DO NOT CALL THIS TWICE FROM INBETWEEN va_start() and va_end()
 * this is not save for old ANSI-Standard!
 *******************************************************************************/
static void write_err(const char *prefix, const char *fmt, va_list arg) {
    if (have_logfile() || check_errstream() != NULL) {
        FILE *f = get_errstream();
        ndm_print(f, PRINT_ALWAYS, prefix);
        ndm_printlist(f, PRINT_ALWAYS, fmt, arg);
        ndm_print(f, PRINT_ALWAYS, "\n\n");
        fflush(f);
    }
} /** write_err() **/

/*******************************************************************************
 *
 *******************************************************************************/
int _warning(const char *fmt, ...) {
    va_list arg;

    va_start(arg, fmt);
    write_out("WARNING:\n", fmt, arg);
    va_end(arg);

    va_start(arg, fmt);
    write_err("WARNING:\n", fmt, arg);
    va_end(arg);

    return 0;
} /** _warning() **/

/*******************************************************************************
 *
 *******************************************************************************/
static void terminate_now(void) {
    ndm_parallel_abort();
    exit(EXIT_FAILURE);
} /** terminate_now() **/

/*******************************************************************************
 *
 *******************************************************************************/
int terminate(const char *fmt, ...) {
    va_list arg;

    va_start(arg, fmt);
    write_out("ERROR:\n", fmt, arg);
    va_end(arg);

    va_start(arg, fmt);
    write_err("ERROR:\n", fmt, arg);
    va_end(arg);

    terminate_now();
    return 0;
} /** terminate() **/

/*******************************************************************************
 *
 *******************************************************************************/
int ndm_exit(const char *fmt, ...) {
    va_list arg;

    va_start(arg, fmt);
    write_out("ERROR:\n", fmt, arg);
    va_end(arg);

    va_start(arg, fmt);
    write_err("ERROR:\n", fmt, arg);
    va_end(arg);

    terminate_now();
    return 0;
} /**  ndm_exit() **/

/*******************************************************************************
 *
 *******************************************************************************/
int schedule_terminate(const char *fmt, ...) {
    va_list arg;

    if (mpi_ndomains() == 1) {
        va_start(arg, fmt);
        write_out("ERROR:\n", fmt, arg);
        va_end(arg);

        va_start(arg, fmt);
        write_err("ERROR:\n", fmt, arg);
        va_end(arg);

        terminate_now();
    } else
        /*----------------------------------------------------------------------------
        | print 1st error message only, suppress others until exit!
        ----------------------------------------------------------------------------*/
        if (error_flag == 0) {
        va_start(arg, fmt);
        write_out("ERROR:\n", fmt, arg);
        va_end(arg);

        va_start(arg, fmt);
        write_err("ERROR:\n", fmt, arg);
        va_end(arg);

        error_flag = 1;
    }

    return 0;
} /** schedule_terminate() **/

/*******************************************************************************
 *
 *******************************************************************************/
int fatal_error_check(void) {
    const int ndom = mpi_ndomains();
    int i = error_flag ? mpi_thisdomain() : ndom;
    myparallel_int_globalmin(&i, 1, 4711);
    if (i < ndom) {
        ndm_print(get_outstream(), PRINT_ALWAYS,
                "\nFATAL ERROR at least for domain (MPI rank) %d\n\n", i);
        if (have_logfile() || check_errstream() != NULL)
            ndm_print(get_errstream(), PRINT_ALWAYS,
                "\nFATAL ERROR at least for domain (MPI rank) %d\n\n", i);
        if (mpi_thisdomain() == i)
            ndm_print(stderr, PRINT_ALWAYS,
                "\nFATAL ERROR at least for domain (MPI rank) %d\n\n", i);
        fflush(stderr);
        mpi_parallel_end();
        ndm_print_close_stdout_stderr();
        exit(EXIT_FAILURE);
    }
    return 0;
}/** fatal_error_check() **/

