#ifndef TIMER_H00
#define TIMER_H00

const int MAP_SIZE = 100;

/**
 * Must be called before any other procedure.
 */
void init_timer();

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
double end_cpu_time(int id);

#endif
