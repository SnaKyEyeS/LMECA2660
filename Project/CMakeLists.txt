# Minimum version
cmake_minimum_required(VERSION 3.10.2)

# Project name
project(Project CXX)

# Exec name
set(EXEC_NAME project)

# CMake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Math
set(CMAKE_REQUIRED_LIBRARIES m)
unset(CMAKE_REQUIRED_LIBRARIES)
set(LIBM_LIBRARIES m)

# PETSc
set(PETSC_DIR ${PROJECT_SOURCE_DIR}/petsc-3.12.4)
set(PETSC_ARCH arch-linux2-c-debug)
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})


# Project files
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/inc)
file(GLOB C_SRCS
    "${PROJECT_SOURCE_DIR}/src/*.c"
	"${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Use G++ for compiling/linking the source files
set_source_files_properties(${C_SRCS} PROPERTIES LANGUAGE CXX)
add_executable(${EXEC_NAME} ${C_SRCS})
target_link_libraries(${EXEC_NAME} ${LIBM_LIBRARIES} ${PETSC_LIBRARIES})
