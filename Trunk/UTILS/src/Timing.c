/*******************************************************************************
 * File:        Timing.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#include "Timing.h"
#include "Utils.h"
#include "Log.h"

/*------------------------------------------------------------------------------
| static variables and structs
------------------------------------------------------------------------------*/
static NDMDouble user_time(void);

static NDMDouble system_time(void);

static Stopwatch watch;

static int output_timing_level = -1;

static AllTimes old_times = {0.0, 0.0, 0.0};

/*******************************************************************************
 *
 ******************************************************************************/
Stopwatch initialise_timing_output(int output_onoff) {
    if (output_timing_level != -1) {
        output_timing_level = output_onoff;
        return watch;
    }

    output_timing_level = output_onoff;

    init_stopwatch(&(watch));

    cont_stopwatch(&(watch));

    return watch;

} /** initialise_timing_output **/

/*******************************************************************************
 *
 ******************************************************************************/
void reset_timing(void) {
    if (output_timing_level < 0)
        initialise_timing_output(1);

    set_current_times(&old_times);
} /** reset_timing() **/

/*******************************************************************************
 *
 ******************************************************************************/
void initialise_timing(void) {
    reset_timing();
} /** initialise_timing() **/

/*******************************************************************************
 *
 ******************************************************************************/
void set_current_times(AllTimes *all_times) {
    if (output_timing_level < 0)
        initialise_timing_output(1);

    all_times->r = read_stopwatch(&(watch));
    all_times->u = user_time();
    all_times->s = system_time();

} /** set_current_times() **/

/*******************************************************************************
 *
 ******************************************************************************/
void output_timing(char *label) {
    output_this_time(label);
    output_total_time();
} /** output_timing() **/

/*******************************************************************************
 *
 ******************************************************************************/
void output_total_time(void) {
    AllTimes total_times;

    if (output_timing_level < 0)
        initialise_timing_output(1);
    else
        if (!output_timing_level)
        return;

    set_current_times(&total_times);

    ndm_msg(" total time in seconds: real %9.3e, user %9.3e sys %9.3e\n",
            total_times.r, total_times.u, total_times.s);

} /** output_total_time() **/

/*******************************************************************************
 *
 ******************************************************************************/
void output_this_time(char *label) {
    AllTimes this_times;
    NDMDouble du, dr, ds;

    if (output_timing_level < 0)
        initialise_timing_output(1);
    else
        if (!output_timing_level)
        return;

    set_current_times(&this_times);

    dr = this_times.r - old_times.r;
    du = this_times.u - old_times.u;
    ds = this_times.s - old_times.s;

    old_times.r = this_times.r;
    old_times.u = this_times.u;
    old_times.s = this_times.s;

    ndm_msg("Timing for: '%s'\n", label);
    ndm_msg("  this time in seconds: real %9.3e, user %9.3e sys %9.3e\n",
            dr, du, ds);

} /** output_this_time() **/

/*******************************************************************************
 * get the user time since last call
 ******************************************************************************/
NDMDouble get_user_time(AllTimes *all_times) {
    NDMDouble d_u_time = user_time() - all_times->u;

    all_times->u += d_u_time;

    return d_u_time;

} /** get_user_time() **/

/*******************************************************************************
 * get the real time since last call
 ******************************************************************************/
NDMDouble get_real_time(AllTimes *all_times) {
    NDMDouble d_r_time = read_stopwatch(&(watch)) - all_times->r;

    all_times->r += d_r_time;

    return d_r_time;

} /** get_real_time() **/

/*******************************************************************************
 * get the system time since last call
 ******************************************************************************/
NDMDouble get_system_time(AllTimes *all_times) {
    NDMDouble d_s_time = system_time() - all_times->s;

    all_times->s += d_s_time;

    return d_s_time;

} /** get_system_time() **/

/*******************************************************************************
 * get elapsed user time
 ******************************************************************************/
static NDMDouble user_time(void) {

#ifndef RUSAGE_SELF
    {
        return 0;
    }
#endif
#ifdef VECTORCOMPUTER
    {
        return 0;
    }
#else
    {
        typedef struct rusage rusage_type;
        rusage_type r;

        if (getrusage(RUSAGE_SELF, &r) != 0)
            ndm_msg("Cannot get system time\n");
        return r.ru_utime.tv_sec + 1e-6 * r.ru_utime.tv_usec;
    }
#endif
} /** user_time() **/

/*******************************************************************************
 * get elapsed system time
 ******************************************************************************/
static NDMDouble system_time(void) {
#ifndef RUSAGE_SELF
    {
        return 0;
    }
#endif
#ifdef VECTORCOMPUTER
    {
        return 0;
    }
#else
    {
        typedef struct rusage rusage_type;
        rusage_type r;

        if (getrusage(RUSAGE_SELF, &r) != 0)
            ndm_msg("Cannot get system time\n");
        return r.ru_stime.tv_sec + 1e-6 * r.ru_stime.tv_usec;
    }
#endif

} /** system_time() **/

