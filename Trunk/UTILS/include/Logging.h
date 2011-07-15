/*******************************************************************************
 * File:        Logging.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _LOGGING_H
#define _LOGGING_H

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*******************************************************************************
* Open logfiles an redirect the error channel and the output channel to files.
* no logfile-redirection for prefix == NULL
* redirect to /dev/null for thisdom > ndom!
*******************************************************************************/
void open_logfiles(const char *prefix, const char *module,
                   int thisdom, int ndom);

/*******************************************************************************
* redirect high process numbers (thisdom >= ndom) to /dev/null
*******************************************************************************/
void stdout2dev0(int thisdom, int ndom);

/*******************************************************************************
* check if logfiles exist
*******************************************************************************/
int have_logfile(void);

/*******************************************************************************
*
*******************************************************************************/
void close_logfiles(void);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _LOGGING_H  */
