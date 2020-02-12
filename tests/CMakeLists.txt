add_executable(testio io.cc)
add_executable(testimpute_r2 impute_r2.cc)
include_directories(".")
include_directories("../lib")


add_test(io ${EXECUTABLE_OUTPUT_PATH}/testio)
add_test(imputer2 ${EXECUTABLE_OUTPUT_PATH}/testimpute_r2)


target_link_libraries(testio geltools_)
target_link_libraries(testimpute_r2 geltools_)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(testio ${Boost_LIBRARIES})
target_link_libraries(testimpute_r2 ${Boost_LIBRARIES})
