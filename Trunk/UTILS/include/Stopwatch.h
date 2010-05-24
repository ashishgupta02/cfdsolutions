/*******************************************************************************
 * File:        Stopwatch.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _STOPWATCH_H
#define _STOPWATCH_H

#include <time.h>

#include "NDM_TypeDefs.h"

/*------------------------------------------------------------------------------
| Keep C++ compilers from getting confused
------------------------------------------------------------------------------*/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| This module provides a timing functionality within the C-code.
| Please do not use the Stopwatch struct components directly outside this module
| because they might change in future.
| Handle the watch via the function calls declared below.
------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
| Stopwatch structure:
|  real_time   : time in seconds the timer was running till the last stop-call
|  time_start  : absolute time of the last continue-call
|  user_time   : time in seconds the timer was running till the last stop-call
|  clock_start : clock count of the last continue-call
|  running     : binary flag indicating the running status
------------------------------------------------------------------------------*/
typedef struct
{
  NDMDouble real_time;
  NDMDouble time_start;
  NDMDouble user_time;
  clock_t   clock_start;
  int       running;
} Stopwatch;

/*******************************************************************************
* Initializes a stopwatch. Times are set to zero.
* After this call the timer is still stopped.
*******************************************************************************/
extern void init_stopwatch(Stopwatch *stopwatch);

/*******************************************************************************
* Starts/continues a stopped stopwatch but has no effect on a running stopwatch.
*******************************************************************************/
extern void cont_stopwatch(Stopwatch *stopwatch);

/*------------------------------------------------------------------------------
| To combine initialize and continue:
------------------------------------------------------------------------------*/
#define start_stopwatch(sw) (init_stopwatch(sw), cont_stopwatch(sw))


/*******************************************************************************
* Stops a running stopwatch but has no effect on an already stopped stopwatch.
*******************************************************************************/
extern void stop_stopwatch(Stopwatch *stopwatch);

/*******************************************************************************
* Returns the real (or wall clock) time in seconds the stopwatch is/was running.
* The timer is not stopped or continued by this function.
*******************************************************************************/
extern NDMDouble read_stopwatch(Stopwatch *stopwatch);

/*******************************************************************************
* Returns the user (or cpu) time in seconds the stopwatch is/was running.
* The timer is not stopped or continued by this function.
*******************************************************************************/
extern NDMDouble read_stopwatch_user(Stopwatch *stopwatch);

/*******************************************************************************
* msleep(x) is supposed to return after x milliseconds. If appropriate timing
* is not supported on a system, however, msleep() may return immediately.
*******************************************************************************/
void msleep(double ms);

/*------------------------------------------------------------------------------
| Keep C++ compilers from getting confused
------------------------------------------------------------------------------*/
#if defined __cplusplus
}
#endif

#endif /* _STOPWATCH_H */
