/*******************************************************************************
 * File:        Command_Line.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _COMMAND_LINE_H
#define _COMMAND_LINE_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*******************************************************************************
*
*******************************************************************************/
char *ndm_commandline(int argc, char **argv,
                      char **logprefix, char **paramfilename,
                      char *module_name, int dom, int ndom);

char *get_commandline(int argc, char **argv,
                      char **logprefix, char **paramfilename);

char *get_and_print_commandline(int argc, char **argv);

const char *return_commandline(void);
void free_commandline(void);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _COMMAND_LINE_H  */

