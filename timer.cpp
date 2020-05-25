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
	for (int i = 0; i < MAP_SIZE; ++i)
	{
		id_map[i] = -1;
	}
}

/**
 * Start the timer. The id must be between 0 and MAP_SIZE to be valid.
 * Negative ids are ignored.
 * Returns the id that was used which can be passed to end_cpu_time to
 * get elapsed time.
 */
int start(int id)
{
	int size = cpu_clock_starts.size();
	if (id >= MAP_SIZE)
	{
		throw std::invalid_argument(std::string("invalid id, must be between 0 and ") 
				+ std::to_string(MAP_SIZE));
	}

	// get a valid id
	if (id < 0)
	{
		for (int i = 0; i < MAP_SIZE; ++i)
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

	std::clock_t cstart = std::clock();
	int index = id_map[id];

	if (id_map[id] < 0 || index >= size)
	{
		cpu_clock_starts.push_back(cstart);
		id_map[id] = size;
	}
	else
	{
		cpu_clock_starts[index] = cstart;
	}
	return id;
}

/**
 * Returns elapsed time in ms.
 * The id is the integer returned by start or the same one that was passed to it.
 */
double end_cpu_time(int id)
{
	std::clock_t cend = std::clock();

	int size = cpu_clock_starts.size();
	if (id >= MAP_SIZE)
	{
		throw std::invalid_argument(std::string("invalid id, must be between 0 and ") 
				+ std::to_string(MAP_SIZE));
	}

	int index = id_map[id];
	if (index < 0 || index >= size)
	{
		std::cout << "boo " << index << " " << id;
		return 0.0;
	}

	// in POSIX CLOCKS_PER_SEC = 1000000
	double elapsed_time = (cend-cpu_clock_starts[index])/1000.0;
	return elapsed_time;
}
