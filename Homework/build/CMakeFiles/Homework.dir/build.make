# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build

# Include any dependencies generated for this target.
include CMakeFiles/Homework.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Homework.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Homework.dir/flags.make

CMakeFiles/Homework.dir/src/solver/solver.c.o: CMakeFiles/Homework.dir/flags.make
CMakeFiles/Homework.dir/src/solver/solver.c.o: ../src/solver/solver.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/Homework.dir/src/solver/solver.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/Homework.dir/src/solver/solver.c.o   -c /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/solver/solver.c

CMakeFiles/Homework.dir/src/solver/solver.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Homework.dir/src/solver/solver.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/solver/solver.c > CMakeFiles/Homework.dir/src/solver/solver.c.i

CMakeFiles/Homework.dir/src/solver/solver.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Homework.dir/src/solver/solver.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/solver/solver.c -o CMakeFiles/Homework.dir/src/solver/solver.c.s

CMakeFiles/Homework.dir/src/solver/solver.c.o.requires:

.PHONY : CMakeFiles/Homework.dir/src/solver/solver.c.o.requires

CMakeFiles/Homework.dir/src/solver/solver.c.o.provides: CMakeFiles/Homework.dir/src/solver/solver.c.o.requires
	$(MAKE) -f CMakeFiles/Homework.dir/build.make CMakeFiles/Homework.dir/src/solver/solver.c.o.provides.build
.PHONY : CMakeFiles/Homework.dir/src/solver/solver.c.o.provides

CMakeFiles/Homework.dir/src/solver/solver.c.o.provides.build: CMakeFiles/Homework.dir/src/solver/solver.c.o


CMakeFiles/Homework.dir/src/solver/thomas.c.o: CMakeFiles/Homework.dir/flags.make
CMakeFiles/Homework.dir/src/solver/thomas.c.o: ../src/solver/thomas.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/Homework.dir/src/solver/thomas.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/Homework.dir/src/solver/thomas.c.o   -c /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/solver/thomas.c

CMakeFiles/Homework.dir/src/solver/thomas.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Homework.dir/src/solver/thomas.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/solver/thomas.c > CMakeFiles/Homework.dir/src/solver/thomas.c.i

CMakeFiles/Homework.dir/src/solver/thomas.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Homework.dir/src/solver/thomas.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/solver/thomas.c -o CMakeFiles/Homework.dir/src/solver/thomas.c.s

CMakeFiles/Homework.dir/src/solver/thomas.c.o.requires:

.PHONY : CMakeFiles/Homework.dir/src/solver/thomas.c.o.requires

CMakeFiles/Homework.dir/src/solver/thomas.c.o.provides: CMakeFiles/Homework.dir/src/solver/thomas.c.o.requires
	$(MAKE) -f CMakeFiles/Homework.dir/build.make CMakeFiles/Homework.dir/src/solver/thomas.c.o.provides.build
.PHONY : CMakeFiles/Homework.dir/src/solver/thomas.c.o.provides

CMakeFiles/Homework.dir/src/solver/thomas.c.o.provides.build: CMakeFiles/Homework.dir/src/solver/thomas.c.o


CMakeFiles/Homework.dir/src/main.c.o: CMakeFiles/Homework.dir/flags.make
CMakeFiles/Homework.dir/src/main.c.o: ../src/main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/Homework.dir/src/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/Homework.dir/src/main.c.o   -c /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/main.c

CMakeFiles/Homework.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Homework.dir/src/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/main.c > CMakeFiles/Homework.dir/src/main.c.i

CMakeFiles/Homework.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Homework.dir/src/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/src/main.c -o CMakeFiles/Homework.dir/src/main.c.s

CMakeFiles/Homework.dir/src/main.c.o.requires:

.PHONY : CMakeFiles/Homework.dir/src/main.c.o.requires

CMakeFiles/Homework.dir/src/main.c.o.provides: CMakeFiles/Homework.dir/src/main.c.o.requires
	$(MAKE) -f CMakeFiles/Homework.dir/build.make CMakeFiles/Homework.dir/src/main.c.o.provides.build
.PHONY : CMakeFiles/Homework.dir/src/main.c.o.provides

CMakeFiles/Homework.dir/src/main.c.o.provides.build: CMakeFiles/Homework.dir/src/main.c.o


# Object files for target Homework
Homework_OBJECTS = \
"CMakeFiles/Homework.dir/src/solver/solver.c.o" \
"CMakeFiles/Homework.dir/src/solver/thomas.c.o" \
"CMakeFiles/Homework.dir/src/main.c.o"

# External object files for target Homework
Homework_EXTERNAL_OBJECTS =

Homework: CMakeFiles/Homework.dir/src/solver/solver.c.o
Homework: CMakeFiles/Homework.dir/src/solver/thomas.c.o
Homework: CMakeFiles/Homework.dir/src/main.c.o
Homework: CMakeFiles/Homework.dir/build.make
Homework: CMakeFiles/Homework.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable Homework"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Homework.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Homework.dir/build: Homework

.PHONY : CMakeFiles/Homework.dir/build

CMakeFiles/Homework.dir/requires: CMakeFiles/Homework.dir/src/solver/solver.c.o.requires
CMakeFiles/Homework.dir/requires: CMakeFiles/Homework.dir/src/solver/thomas.c.o.requires
CMakeFiles/Homework.dir/requires: CMakeFiles/Homework.dir/src/main.c.o.requires

.PHONY : CMakeFiles/Homework.dir/requires

CMakeFiles/Homework.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Homework.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Homework.dir/clean

CMakeFiles/Homework.dir/depend:
	cd /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build /home/jyl/Documents/GitHub/SnaKyEyeS.github.io/LMECA2660/Homework/build/CMakeFiles/Homework.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Homework.dir/depend

