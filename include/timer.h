/*
 * Contains functions for easily timing sections of code.
 * Usage:
 *   1. Invoke init_timer once for initializaiton
 *
 *   2. Invoke start and keep track of id returned before
 *      the section of code that needs to be timed
 *
 *   3. At the end, invoking end_cpu_time will return
 *      approximately how much cpu time elapsed since
 *      start was called
 */

#ifndef TIMER_H00
#define TIMER_H00
#include <vector>

const size_t MAP_SIZE = 100;

struct time_metric 
{
    const double cpu_time;
    const double wall_time;
};

/**
 */
void init_timer();

/**
 * Start the timer. 
 * Returns the id that was used which can be passed to end_cpu_time to
 * get elapsed time.
 */
int start();

/**
 * Start the timer. The id must be between 0 and MAP_SIZE to be valid.
 * Negative ids are ignored.
 * Returns the id that was used which can be passed to end_cpu_time to
 * get elapsed time.
 */
int start(int id);

/**
 * Returns elapsed time in ms.
 * The id is the integer returned by start or the same one that was passed to it.
 */
const time_metric end(int id);

/**
 * Prints out the cpu and wall time to stdout;
 */
void printTime(const time_metric &t);

std::ostream& operator<<(std::ostream &os, const time_metric &t);

#endif
