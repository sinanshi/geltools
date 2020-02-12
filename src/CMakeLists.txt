file(GLOB src_files "*.cpp")
file(GLOB header_files "*.h")
include_directories(".")
include_directories("../lib")

add_executable(geltools  ${src_files})
target_link_libraries(geltools ${Boost_LIBRARIES}) #${BGEN_LIBRARIES})
target_link_libraries(geltools ${HTSLIB_LIBRARY}) #${BGEN_LIBRARIES})


add_library(geltools_ SHARED ${src_files})
target_link_libraries(geltools_ ${Boost_LIBRARIES}) #${BGEN_LIBRARIES})
target_link_libraries(geltools_ ${HTSLIB_LIBRARY}) #${BGEN_LIBRARIES})
