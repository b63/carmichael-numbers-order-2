link_libraries(pthread m gmp ntl)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/strategy_2")

add_executable(gen_nonrigid ../counting_factors.cpp ../timer.cpp
               ../util.cpp nonrigid.cpp gen_nonrigid.cpp
               )

add_executable(all_nonrigid_pairs ../timer.cpp ../util.cpp 
              nonrigid.cpp all_possible_nonrigid_pairs.cpp
              )
