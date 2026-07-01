/*******************************************************************************
 * File:        Stopwatch.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "Stopwatch.h"
#include "Utils.h"

/*****************************************************************************
 * system_time() returns the wallclock time in seconds (relative to some fixed
 * point in the past, e.g. in POSIX 1.1.1970, 00:00:00 UTC)
 ****************************************************************************/
static double system_time(void) {
#if defined(_WIN32)
    LARGE_INTEGER frequency;
    LARGE_INTEGER counter;
    if (!QueryPerformanceFrequency(&frequency) || !QueryPerformanceCounter(&counter))
        return 0.0;
    return (double)counter.QuadPart / (double)frequency.QuadPart;
#else
    struct timeval tv;
    if (gettimeofday(&tv, NULL) != 0)
        return 0.0;
    else
        return (double) tv.tv_sec + 1e-6 * (double) tv.tv_usec;
#endif
}

/*******************************************************************************
 *
 *******************************************************************************/
void msleep(double ms) {
    double start = system_time();
    if (start > 0.001) while (system_time() < start + ms * 0.001) {
        ;
    }
}

/*******************************************************************************
 *
 *******************************************************************************/
void init_stopwatch(Stopwatch *sw) {
    sw->user_time = 0.0;
    sw->real_time = 0.0;
    sw->running = FALSE;

} /** init_stopwatch() **/

/*******************************************************************************
 *
 *******************************************************************************/
void cont_stopwatch(Stopwatch *sw) {
    if (sw->running)
        return;

    sw->time_start = system_time();
    sw->clock_start = clock();
    sw->running = TRUE;

} /** cont_stopwatch() **/

/*******************************************************************************
 *
 *******************************************************************************/
void stop_stopwatch(Stopwatch *sw) {
    if (!sw->running)
        return;

    sw->real_time += system_time() - sw->time_start;
    sw->user_time += (clock() - sw->clock_start) / (CLOCKS_PER_SEC * 1.);
    sw->running = FALSE;

} /** stop_stopwatch() **/

/*******************************************************************************
 *
 *******************************************************************************/
NDMDouble read_stopwatch(Stopwatch *sw) {
    if (!sw->running)
        return sw->real_time;

    return sw->real_time + system_time() - sw->time_start;

} /** read_stopwatch() **/

/*******************************************************************************
 *
 *******************************************************************************/
NDMDouble read_stopwatch_user(Stopwatch *sw) {
    if (!sw->running)
        return sw->user_time;

    return sw->user_time + (clock() - sw->clock_start) / (CLOCKS_PER_SEC * 1.);

} /** read_stopwatch_user() **/

