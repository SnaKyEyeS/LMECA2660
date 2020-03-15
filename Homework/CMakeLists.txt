# Minimum version
cmake_minimum_required(VERSION 3.7)

# Project name
project(Homework)

# Exec name
set(EXEC_NAME Homework)

# CMake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Math
set(CMAKE_REQUIRED_LIBRARIES m)
unset(CMAKE_REQUIRED_LIBRARIES)
set(LIBM_LIBRARIES m)

# Project files
include_directories(${PROJECT_SOURCE_DIR}/src/solver ${PROJECT_SOURCE_DIR}/src)
file(GLOB SRCS
    "${PROJECT_SOURCE_DIR}/src/solver/*.h"
	"${PROJECT_SOURCE_DIR}/src/solver/*.c"
	"${PROJECT_SOURCE_DIR}/src/*.h"
	"${PROJECT_SOURCE_DIR}/src/*.c"
)

add_executable(${EXEC_NAME} ${SRCS})
target_link_libraries(${EXEC_NAME} ${LIBM_LIBRARIES})
