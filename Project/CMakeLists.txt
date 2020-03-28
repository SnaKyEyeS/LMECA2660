# Minimum version
cmake_minimum_required(VERSION 3.10.2)

# Project name
project(Project)

# Exec name
set(EXEC_NAME Project)

# CMake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Math
set(CMAKE_REQUIRED_LIBRARIES m)
unset(CMAKE_REQUIRED_LIBRARIES)
set(LIBM_LIBRARIES m)

# PETSc
find_package(PETSc 3.12.4 REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

# Project files
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/inc)
file(GLOB SRCS
    "${PROJECT_SOURCE_DIR}/inc/*.h"
	"${PROJECT_SOURCE_DIR}/inc/*.c"
	"${PROJECT_SOURCE_DIR}/src/*.h"
	"${PROJECT_SOURCE_DIR}/src/*.c"
)

add_executable(${EXEC_NAME} ${SRCS})
target_link_libraries(${EXEC_NAME} ${LIBM_LIBRARIES} ${PETSC_LIBRARIES})
