/*
 * Contains functions for easily timing sections of code.
 * Usage:
 *   1. Invoke init_timer once for initializaiton
 *
 *   2. Invoke start and keep track of id returned before
 *      the section of code that needs to be timed
 *
 *   3. At the end, invoking end(id) will return
 *      approximately how much cpu time elapsed since
 *      start was called. get_cpu_time(id) and get_wall_time(id)
 *      can also be called after end(id).
 */
#include <time.h>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <vector>
#include <timer.h>

typedef std::chrono::high_resolution_clock::time_point time_point;
typedef std::chrono::duration<double, std::milli> double_duration;

std::vector<std::clock_t> cpu_clock_starts;
std::vector<time_point> wall_clock_starts;
int id_map[MAP_SIZE];


/**
 * Initialize the timer. Must be called once before any
 * other procedure in this module..
 */
void init_timer()
{
    for (size_t i = 0; i < MAP_SIZE; i++)
    {
        id_map[i] = -1;
    }
}


/**
 * Start the timer. 
 * Returns the id that was used which can be passed to end_cpu_time to
 * get elapsed time.
 */
int start()
{
    return start(-1);
}

/**
 * Start the timer. The id must be between 0 and MAP_SIZE to be valid.
 * Negative ids are ignored.
 * Returns the id that was used which can be passed to end_cpu_time to
 * get elapsed time.
 */
int start(int id)
{
    // get a valid id if id is negative
    if (id < 0)
    {
        for (size_t i = 0; i < MAP_SIZE; i++)
        {
            if (id_map[i] < 0)
            {
                id = i;
                break;
            }

        }
    }

    if (id < 0)
    {
        throw std::invalid_argument(std::string("could not find valid id")) ;
    }

    size_t pos_id { (size_t) id };
    size_t size { cpu_clock_starts.size() };
    if (pos_id >= MAP_SIZE)
    {
        throw std::invalid_argument(std::string("invalid id, must be between 0 and ") 
                + std::to_string(MAP_SIZE));
    }

    time_point wstart 
        { std::chrono::high_resolution_clock::now() };
    std::clock_t cstart { std::clock() };
    int index { id_map[pos_id] };

    if (index < 0 || (size_t) index >= size)
    {
        // negative index means id has not been used
        // index > size means id doesn't correspond to a valid index
        // so assign a new index to id
        cpu_clock_starts.push_back(cstart);
        wall_clock_starts.push_back(wstart);
        id_map[pos_id] = size;
    }
    else
    {
        // valid index, so just overwrite
        cpu_clock_starts[index] = cstart;
        wall_clock_starts[index] = wstart;
    }

    // return the id that was used
    return pos_id;
}

/**
 * Returns elapsed time in seconds as time_metric, which has two members:
 * cpu_time and wall_time. 
 * The id is the integer returned by start or the same one that was passed to it.
 */
const time_metric end(int id)
{
    std::clock_t cend { std::clock() };
    time_point wend { std::chrono::high_resolution_clock::now() };

    size_t size = cpu_clock_starts.size();
    if ((size_t) id >= MAP_SIZE)
    {
        throw std::invalid_argument(std::string("invalid id, must be between 0 and ") 
                + std::to_string(MAP_SIZE));
    }

    size_t index = id_map[id];
    if (index >= size)
    {
        // invalid index
        return time_metric {-1, -1};
    }

    // in POSIX CLOCKS_PER_SEC = 1 000 000
    double elapsed_time_cpu { (cend-cpu_clock_starts[index])/1000.0 };
    double elapsed_time_wall { std::chrono
        ::duration_cast<double_duration>(wend - wall_clock_starts[index]).count() };

    return time_metric {elapsed_time_cpu, elapsed_time_wall};
}


std::ostream&
operator<<(std::ostream &os, const time_metric &t)
{
    os << "(" << t.wall_time << ")";
    return os;
}


void printTime(const time_metric &t)
{
    std::cout << "(cpu time " << t.cpu_time << "ms, wall time " 
        <<  t.wall_time << "ms)\n";
}
