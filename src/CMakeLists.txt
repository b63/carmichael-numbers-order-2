link_libraries(pthread m gmp ntl)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(generate_cprimes generate_cprimes.cpp util.cpp)

add_executable(generate_cprimes_order_2 
              generate_cprimes_order_2.cpp timer.cpp util.cpp
              counting_factors.cpp primality.cpp construct_P.cpp
              )

add_executable(construct_P construct_P_main.cpp construct_P.cpp
              timer.cpp util.cpp counting_factors.cpp primality.cpp
              )

add_subdirectory(strategy_2)
