/*******************************************************************************
 * File:        Command_Line.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Check_Malloc.h"
#include "Command_Line.h"
#include "Utils.h"
#include "Log.h"
#include "Warn_Error.h"
#include "Logging.h"

static char *commandline = NULL;

/*******************************************************************************
 *
 ******************************************************************************/
const char *return_commandline(void) {
    const char *my_commandline = commandline != NULL ? commandline : " ";

    return my_commandline;
} /** return_commandline() **/

/*******************************************************************************
 *
 ******************************************************************************/
void free_commandline(void) {
    check_free(commandline);
} /** free_commandline() **/

/*******************************************************************************
 *
 ******************************************************************************/
char *ndm_commandline(int argc, char **argv,
        char **logprefix, char **paramfilename,
        char *module_name, int dom, int ndom) {
    char *cline;

    /*--------------------------------------------------------------------------
    | commandline
    --------------------------------------------------------------------------*/
    cline = get_commandline(argc, argv, logprefix, paramfilename);

    /*--------------------------------------------------------------------------
    | redirect stdout/stderr
    --------------------------------------------------------------------------*/
    open_logfiles(*logprefix, module_name, dom, ndom);

    return cline;
} /** ndm_commandline() **/

/*******************************************************************************
 *
 ******************************************************************************/
char *get_commandline(int argc, char **argv,
        char **logprefix, char **paramfilename) {
    *logprefix = NULL;

    if (argc < 2)
        terminate("Usage: %s <paramfile> <logfile_prefix>", argv[0]);

    /*--------------------------------------------------------------------------
    | assign arguments
    --------------------------------------------------------------------------*/
    *paramfilename = argv[1];
    if (argc >= 3)
        *logprefix = argv[2];

    /*--------------------------------------------------------------------------
    | commandline
    --------------------------------------------------------------------------*/
    commandline = get_and_print_commandline(argc, argv);
    return commandline;
} /** get_commandline() **/

/*******************************************************************************
 * store and print commandline
 ******************************************************************************/
char *get_and_print_commandline(int argc, char **argv) {
    int i, len;
    char *cline;

    /*--------------------------------------------------------------------------
    | print commandline
    --------------------------------------------------------------------------*/
    for (i = 0, len = 1; i < argc; i++)
        len += (int) strlen(argv[i]) + 1;

    cline = (char *) check_malloc(len * sizeof (char));

    for (i = 0, len = 0; i < argc; i++) {
        sprintf(&cline[len], "%s%c", argv[i], (i == argc - 1) ? '\n' : ' ');
        len += (int) strlen(argv[i]) + 1;
    }
    ndm_msg("\nCommandline:\n%s\n", cline);

    return cline;
} /** get_and_print_commandline() **/

