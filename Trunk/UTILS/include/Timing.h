/*******************************************************************************
 * File:        Timing.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _TIMING_H
#define _TIMING_H

#include "Stopwatch.h"

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| container for all times
------------------------------------------------------------------------------*/
typedef struct AllTimes_struct {
    NDMDouble r; /* real time    */
    NDMDouble u; /* user time    */
    NDMDouble s; /* system time  */
} AllTimes;

/*******************************************************************************
* initial call with on/off switch for output
*******************************************************************************/
Stopwatch initialise_timing_output(int output_onoff);

/*******************************************************************************
* call to init timer (output is switched on)
*******************************************************************************/
void initialise_timing(void);

void reset_timing(void);

/*******************************************************************************
* calls to get times
*******************************************************************************/
void current_times(NDMDouble *real_t, NDMDouble *user_t, NDMDouble *sys_t);

/*******************************************************************************
* calls to output timing
*******************************************************************************/
void output_timing(char *label);

void output_total_time(void);

void output_this_time(char *label);

/*******************************************************************************
*
*******************************************************************************/
void set_current_times(AllTimes *all_times);

NDMDouble get_user_time(AllTimes *all_times);

NDMDouble get_system_time(AllTimes *all_times);

NDMDouble get_real_time(AllTimes *all_times);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _TIMING_H */

