cmake_minimum_required(VERSION 3.8)

project(carmichael)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)

set(GLOBAL_LINK_FLAGS "-lpthread -lm -lgmp -lntl")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLOBAL_CXX_FLAGS}")

include_directories(BEFORE ${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
