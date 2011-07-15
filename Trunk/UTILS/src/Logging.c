/*******************************************************************************
 * File:        Logging.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <string.h>

#include "Utils.h"
#include "Log.h"
#include "Warn_Error.h"
#include "Logging.h"

#define DEVNULL "/dev/null"

static int have_stderr = 0;

static void get_name(char *name,
        const char *prefix, const char *module, char *stderr_out,
        int thisdom, int ndom);

/*******************************************************************************
 * New method for loging in stdour/stderr-files without general redirect using
 * log-functions
 ******************************************************************************/
static void open_logfiles_new(const char *prefix, const char *module,
        int thisdom, int ndom) {
    if (prefix != NULL) {
        char filename[MAX_FILENAME_LEN];

        get_name(filename, prefix, module, "stdout", thisdom, ndom);
        ndm_print_set_stdout_filename(filename);
        get_name(filename, prefix, module, "stderr", thisdom, ndom);
        ndm_print_set_stderr_filename(filename);
    }
} /** open_logfiles_new() **/

/*******************************************************************************
 * OLD METHOD should not longer used
 * Open logfiles an redirect the error channel and the output channel to files.
 * redirect to /dev/null for thisdom > ndom!
 *
 * ToDo
 * remove this function and have_logfile()-function if new method is established
 ******************************************************************************/
static void open_logfiles_old(const char *prefix, const char *module,
        int thisdom, int ndom) {
    char filename[MAX_FILENAME_LEN];
    strcpy(filename, DEVNULL);

    if (thisdom < ndom) {
        if (prefix != NULL) {
            get_name(filename, prefix, module, "stdout", thisdom, ndom);
            if (freopen(filename, "w", stdout) == NULL)
                terminate("Cannot redirect 'stdout' to file '%s'!", filename);


            get_name(filename, prefix, module, "stderr", thisdom, ndom);
            if (freopen(filename, "w", stderr) == NULL)
                terminate("Cannot redirect 'stderr' to file '%s'!", filename);

            have_stderr = 1;
        }
    } else {
        stdout2dev0(thisdom, ndom);
    }

} /** open_logfiles_old() **/

/*******************************************************************************
 * redirect stderr & stdout
 ******************************************************************************/
void open_logfiles(const char *prefix, const char *module,
        int thisdom, int ndom) {
    int use_new = TRUE;
    if (use_new)
        open_logfiles_new(prefix, module, thisdom, ndom);
    else
        open_logfiles_old(prefix, module, thisdom, ndom);
} /** open_logfiles() **/

/*******************************************************************************
 * redirect high process numbers (thisdom >= ndom) to /dev/null
 ******************************************************************************/
void stdout2dev0(int thisdom, int ndom) {
    if (thisdom >= ndom) {
        char filename[MAX_FILENAME_LEN];

        strcpy(filename, DEVNULL);
        if (freopen(filename, "w", stdout) == NULL)
            terminate("Cannot redirect 'stdout' to file '%s'!", filename);

        strcpy(filename, DEVNULL);
        if (freopen(filename, "w", stderr) == NULL)
            terminate("Cannot redirect 'stderr' to file '%s'!", filename);
    }
} /** stdout2dev0() **/

/*******************************************************************************
 *
 ******************************************************************************/
int have_logfile(void) {
    return have_stderr;
} /** have_logfile() **/

/*******************************************************************************
 *
 ******************************************************************************/
void close_logfiles(void) {
    ndm_print_close_stdout_stderr(); /* needed if log is in use */
} /** close_logfiles() **/

/*******************************************************************************
 *
 ******************************************************************************/
static void get_name(char *name,
        const char *prefix, const char *module, char *stderr_out,
        int thisdom, int ndom) {
    char suffix[MAX_FILENAME_LEN];
    int slen = strlen(stderr_out) + strlen(module) + 3;
    int sdom = 7;
    int len = (ndom > 1) ? (slen + sdom + 1) : slen;

    if (prefix == NULL)
        fatal_error("got no prefix");
    if (stderr_out == NULL)
        fatal_error("got wrong name");

    if (len >= MAX_FILENAME_LEN)
        fatal_error("logfile-name suffix is too long\n");

    if (ndom > 1)
        sprintf(suffix, ".%s.%s.%d", module, stderr_out, thisdom);
    else
        sprintf(suffix, ".%s.%s", module, stderr_out);

    if (len + strlen(prefix) >= MAX_FILENAME_LEN)
        fatal_error("logfile-name '%s%s' is too long\n", prefix, suffix);

    sprintf(name, "%s%s", prefix, suffix);

    return;
} /** get_name() **/

