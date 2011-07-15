/*******************************************************************************
 * File:        Log.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _LOG_H
#define	_LOG_H

#ifdef	__cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| prints a line (unique design of screen output)
------------------------------------------------------------------------------*/
#define printline ndm_msg(" ---------------------------------------"\
                         "---------------------------------------\n");
#define printshortline ndm_msg("-----------------------------------"\
                              "------------------\n");

/*------------------------------------------------------------------------------
| print function to redirect stdout/stderr - prints
| with checks if domain numbers and output level is active for printing
------------------------------------------------------------------------------*/
# define printf ndm_msg
FILE *check_errstream(void);/*get pointer to stderr or NULL if no logfile set*/
FILE *get_errstream(void); /* get setting for stderr */
FILE *get_outstream(void); /* get setting for stdout */

int ndm_msg(const char *fmt, ...);                /* interface to ndm_print.. */
int ndm_errmsg(const char *fmt, ...);             /* interface to ndm_print.. */

int ndm_print_output(int io_level);           /* return if io_level is active */
int ndm_print(FILE *file, int output, const char *fmt, ...); /*replace printf */
int ndm_print2(int output, const char *fmt, ...);            /*replace printf */
int ndm_printlist(FILE *file, int output, const char *fmt, va_list arg);

void ndm_print_set_output_domain(int max_output_dom);
void ndm_print_set_output_level(int my_output_level);
void ndm_print_set_flushtime(double flushtime);   /* set flushtime in seconds */

/*------------------------------------------------------------------------------
| redirect logoutput
| to change settings (from comandline or own previous settings)
| you have to use close-function 1st before ..set_std-functions!
------------------------------------------------------------------------------*/
void ndm_print_set_stderr_filename(const char *filename);
void ndm_print_set_stdout_filename(const char *filename);
void ndm_print_close_stdout_stderr(void);

/*------------------------------------------------------------------------------
| Debug-Messages
------------------------------------------------------------------------------*/
int have_debug_level(void);
void init_debug_messages(int set_io_level, int set_debug_level);
void debug_print_line_proc(const char *file, int line);

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 100
#endif

#ifndef NO_MESSAGES /*********************************************************/

#define DEBUG_PRINT debug_print_line_proc( __FILE__, __LINE__);

void debug_message(char *entry, char *name, int io_level,
                   const char *file, int line);

#define ENTRY_MESSAGE(a, b) \
  debug_message("ENTER ROUTINE", a, b, __FILE__, __LINE__);
#define EXIT_MESSAGE(a, b)  \
  debug_message("EXIT ROUTINE",  a, b, __FILE__, __LINE__);
#define STATUS_MESSAGE(a, b) \
  debug_message("INSIDE ROUTINE",a, b, __FILE__, __LINE__);

#else /** ifndef NO_MESSAGES: saves about 1% of CPU! *************************/

#define DEBUG_PRINT
#define ENTRY_MESSAGE(a,b)
#define EXIT_MESSAGE(a,b)
#define STATUS_MESSAGE(a,b)

#endif /** ifndef NO_MESSAGES ************************************************/



#ifdef	__cplusplus
}
#endif

#endif	/* _LOG_H */

