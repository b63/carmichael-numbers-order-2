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
#include <time.h>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <vector>
#include "timer.h"

std::vector<std::clock_t> cpu_clock_starts;
int id_map[MAP_SIZE];

void init_timer()
{
	for (size_t i = 0; i < MAP_SIZE; ++i)
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
		for (size_t i = 0; i < MAP_SIZE; ++i)
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

	size_t pos_id = (size_t) id;
	size_t size = cpu_clock_starts.size();
	if (pos_id >= MAP_SIZE)
	{
		throw std::invalid_argument(std::string("invalid id, must be between 0 and ") 
				+ std::to_string(MAP_SIZE));
	}

	std::clock_t cstart = std::clock();
	int index = id_map[pos_id];

	if (index < 0 || (size_t) index >= size)
	{
		// negative index means id has not been used
		// index > size means id doesn't correspond to a valid index
		// so assing a new index to id
		cpu_clock_starts.push_back(cstart);
		id_map[pos_id] = size;
	}
	else
	{
		// valid index, so just overwrite
		cpu_clock_starts[index] = cstart;
	}

	// return the id that was used
	return pos_id;
}

/**
 * Returns elapsed time in ms.
 * The id is the integer returned by start or the same one that was passed to it.
 */
double end_cpu_time(int id)
{
	std::clock_t cend = std::clock();

	size_t size = cpu_clock_starts.size();
	if ((size_t) id >= MAP_SIZE)
	{
		throw std::invalid_argument(std::string("invalid id, must be between 0 and ") 
				+ std::to_string(MAP_SIZE));
	}

	size_t index = id_map[id];
	if (index >= size)
	{
		std::cout << "boo " << index << " " << id;
		return 0.0;
	}

	// in POSIX CLOCKS_PER_SEC = 1000000
	double elapsed_time = (cend-cpu_clock_starts[index])/1000.0;
	return elapsed_time;
}
