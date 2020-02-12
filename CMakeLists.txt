cmake_minimum_required(VERSION 3.2)

project(geltools)
# Mandatory Requirements

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )
find_package(Boost COMPONENTS program_options iostreams filesystem system REQUIRED)
if (Boost_FOUND)
  include_directories(${BOOST_INCLUDE_DIRS})
endif()



find_package(HTSlib REQUIRED)
#include( UseCython )

# With CMake, a clean separation can be made between the source tree and the
# build tree.  When all source is compiled, as with pure C/C++, the source is
# no-longer needed in the build tree.  However, with pure *.py source, the
# source is processed directly.  To handle this, we reproduce the availability
# of the source files in the build tree.
#add_custom_target( ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
#  ${CMAKE_CURRENT_BINARY_DIR}
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )


#set (CMAKE_CXX_FLAGS "-std=c++11 $")
set(CMAKE_CXX_STANDARD 14) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

set(CMAKE_CXX_FLAGS "-msse -msse2 -msse3 -Wall -ftree-vectorize -DAVX -mavx -O3 -fPIC")
enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
#add_subdirectory(doc)

set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)


#set(doxy_main_page README.md)
