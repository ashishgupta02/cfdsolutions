/*******************************************************************************
 * File:        Log.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "Utils.h"
#include "Log.h"
#include "Stopwatch.h"
#include "Output_Level.h"
#include "Interface_MPI.h"

/*******************************************************************************
 * handle all log-output
 *
 * the functions below are to be used instead of printf or fprintf for all
 * stdout or stderr output
 * this allows to redirect stdout/stderr to files or switch on/off output
 * depending on output-level, domain-number or other criteria for which
 * all relevant settings  are under control inside this module
 ******************************************************************************/
static int output_dom = 99999999;
static int output_level = PRINT_MAXIMUM_LEVEL;

static FILE *ioerr = NULL;
static FILE *ioout = NULL;

/*******************************************************************************
 * additional debug stuff
 ******************************************************************************/
static int debug_level = 100;

/*******************************************************************************
 * replace printf, fprintf, vprintf, vfprintf etc
 * by one single ndm functions ndm_printlist() and provide some other
 * interfaces to this function for easy usage
 *
 * This is to switch stdout/stderr to files on selected processes, i.e.
 * full control of log messages via a single function
 ******************************************************************************/

/*------------------------------------------------------------------------------
| Timer for private fflush function
------------------------------------------------------------------------------*/
static Stopwatch flushtimer;
static int stopwatch_is_initialized = FALSE;
static NDMDouble flushtime_next = 0.0;
static NDMDouble flushtime = -1.0; /* negative value to init default */

/*******************************************************************************
 * interfaces to init ndm_print settings e.g. in param_init or via python
 ******************************************************************************/
void ndm_print_set_flushtime(double flushtime_interval) {
    /*--------------------------------------------------------------------------
    | set default values for private fflush function for flushtime in seconds
    | do not flush all at the same time to help our parallel filesystem
    | take a magic number to flush at different times if possible, this might
    | not help because the procs flushes in different files, but it can not harm
    --------------------------------------------------------------------------*/
    if (flushtime < 0. || flushtime_interval < 0.) {
        int mydom = mpi_thisdomain();
        NDMDouble delta_per_core = 1.3578 * (double) (mydom % 32);

        if (mydom == 0)
            flushtime = 6.0; /*smaller values show a bad effect on CASE */
        else
            flushtime = 60.0 + delta_per_core;

    }
    /*--------------------------------------------------------------------------
    | use parameter settings from outside for values for private fflush function
    --------------------------------------------------------------------------*/
    else {
        flushtime = flushtime_interval;
    }
}

/*******************************************************************************
 * introduce time interval for fflush
 * leaving the exception for fllushing always for debugging: output_level = 99!
 ******************************************************************************/
static void private_fflush(FILE *f) {
    fflush(f);

    if (output_level >= 99) {
        fflush(f);
        return;
    }

    if (!stopwatch_is_initialized) {
        if (flushtime < 0.)
            ndm_print_set_flushtime(-1.0);

        init_stopwatch(&flushtimer);
        stopwatch_is_initialized = TRUE;
        flushtime_next = flushtime;
        cont_stopwatch(&flushtimer);
    }

    if (read_stopwatch_user(&flushtimer) >= flushtime_next) {
        fflush(f);
        flushtime_next += flushtime;
    }
}

/*******************************************************************************
 * get functions for stderr/stdout as defined in init-phase
 ******************************************************************************/
FILE *get_errstream(void) {
    if (ioerr == NULL)
        return stderr;
    else
        return ioerr;
}

FILE *check_errstream(void) {
    return ioerr;
}

FILE *get_outstream(void) {
    if (ioout == NULL)
        return stdout;
    else
        return ioout;
}

/*******************************************************************************
 * criteria for output
 ******************************************************************************/
int ndm_print_output(int io_level) {
    if (io_level == PRINT_ALWAYS) /* do it on all domains! */
        return TRUE;
    else if (io_level <= output_level && mpi_thisdomain() <= output_dom)
        return TRUE;

    return FALSE;
}

/*******************************************************************************
 * helper function for seting stderr/stdout-filenname
 ******************************************************************************/
static FILE *set_stderr_stdout(const char *filename, FILE * oldptr) {
    FILE *myptr = NULL;
    char myfilename[MAX_FILENAME_LEN];

    if (mpi_thisdomain() < output_dom)
        strcpy(myfilename, filename);
    else
        strcpy(myfilename, "/dev/null");

    if (oldptr != NULL)
        ndm_print(oldptr, 0, "file %s already open, close before set\n", filename);
    else
        myptr = fopen(myfilename, "w");
    return myptr;
}

/*******************************************************************************
 * interfaces to init logfile, stderr & stdout
 ******************************************************************************/
void ndm_print_set_stderr_filename(const char *filename) {
    FILE *myptr = set_stderr_stdout(filename, ioerr);

    if (myptr == NULL)
        ndm_print(ioerr, 0, "<%s> could not be opened\n", filename);
    else
        ioerr = myptr;
}

/*******************************************************************************
 * interfaces to init logfile, stderr & stdout
 ******************************************************************************/
void ndm_print_set_stdout_filename(const char *filename) {
    FILE *myptr = set_stderr_stdout(filename, ioout);

    if (myptr == NULL)
        ndm_print(ioerr, 0, "<%s> could not be opened\n", filename);
    else
        ioout = myptr;
}

/*******************************************************************************
 * close output channels
 ******************************************************************************/
void ndm_print_close_stdout_stderr(void) {
    if (ioerr != NULL) {
        fclose(ioerr);
        ioerr = NULL;
    }
    if (ioout != NULL) {
        fclose(ioout);
        ioout = NULL;
    }
}

/*******************************************************************************
 * interfaces to init ndm_print settings done in param_init
 ******************************************************************************/
void ndm_print_set_output_domain(int max_output_dom) {
    output_dom = max_output_dom;
}

/*******************************************************************************
 * interfaces to init ndm_print settings done in param_init
 ******************************************************************************/
void ndm_print_set_output_level(int my_output_level) {
    output_level = my_output_level;
}

/*******************************************************************************
 * interface to ndm_printlist()
 ******************************************************************************/
int ndm_msg(const char *fmt, ...) {
    va_list arg;
    int status;

    va_start(arg, fmt);
    status = ndm_printlist(NULL, DEFAULT_PRINTF_LEVEL, fmt, arg);
    va_end(arg);

    return status;
}

/*******************************************************************************
 * interface to ndm_printlist()
 ******************************************************************************/
int ndm_errmsg(const char *fmt, ...) {
    FILE *ndmerr = get_errstream();
    va_list arg;
    int status;

    va_start(arg, fmt);
    status = ndm_printlist(ndmerr, PRINT_ALWAYS, fmt, arg);
    va_end(arg);

    return status;
}

/*******************************************************************************
 * interface to ndm_printlist()
 ******************************************************************************/
int ndm_print(FILE *file, int output, const char *fmt, ...) {
    int status;

    va_list arg;
    va_start(arg, fmt);
    status = ndm_printlist(file, output, fmt, arg);
    va_end(arg);

    return status;
}

/*******************************************************************************
 * interface to ndm_printlist()
 ******************************************************************************/
int ndm_print2(int output, const char *fmt, ...) {
    int status;

    va_list arg;
    va_start(arg, fmt);
    status = ndm_printlist(NULL, output, fmt, arg);
    va_end(arg);
    ndm_print(NULL, output, "\n");
    return status;
}

/*******************************************************************************
 * replacement of printf or fprintf
 * this should be the only function printing log messages to allow for
 * full control via this method
 ******************************************************************************/
int ndm_printlist(FILE *file, int output, const char *fmt, va_list arg) {
    if (ndm_print_output(output)) {
        FILE *f = (file == NULL) ? get_outstream() : file;

        vfprintf(f, fmt, arg);
        private_fflush(f);
    }
    return 0;
}

/*******************************************************************************
 * additionally some small utilities for debugging
 ******************************************************************************/

/*******************************************************************************
 * debug utilty
 ******************************************************************************/
void init_debug_messages(int set_io_level, int set_debug_level) {
    output_level = set_io_level;
    debug_level = set_debug_level;
}

/*******************************************************************************
 * debug utilty
 ******************************************************************************/
void debug_message(char *entry, char *name, int io_level,
        const char *file, int line) {
    int io = (io_level > 0) ? io_level : output_level;

    if (debug_level <= io)
        ndm_print(NULL, PRINT_ALWAYS, "%s %s: FILE %s LINE %d\n",
            entry, name, file, line);
}

/*******************************************************************************
 * debug utilty
 ******************************************************************************/
void debug_print_line_proc(const char *file, int line) {
    NDMDouble dummy = 1.;
    int process = mpi_thisdomain();

    /*--------------------------------------------------------------------------
    | Dummy communication to syncronize processes in parallel mode
    --------------------------------------------------------------------------*/
    myparallel_globalmax(&dummy, 1, 4711);

    /*--------------------------------------------------------------------------
    | print
    --------------------------------------------------------------------------*/
    ndm_print(NULL, PRINT_ALWAYS, "FILE %s LINE %d PROCESS %d\n",
            file, line, process);

}

/*******************************************************************************
 * debug utilty
 ******************************************************************************/
int have_debug_level(void) {
    if (DEBUG_LEVEL_1 > output_level)
        return 0;
    else
        return 1;
}

