# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/xsarm06/Development/cpp_projects/optimization_methods_x64

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build

# Include any dependencies generated for this target.
include unit_tests/CMakeFiles/fibonacci_test.dir/depend.make

# Include the progress variables for this target.
include unit_tests/CMakeFiles/fibonacci_test.dir/progress.make

# Include the compile flags for this target's objects.
include unit_tests/CMakeFiles/fibonacci_test.dir/flags.make

unit_tests/CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.o: unit_tests/CMakeFiles/fibonacci_test.dir/flags.make
unit_tests/CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.o: ../unit_tests/om_fibonacci_t.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unit_tests/CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.o"
	cd /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests && /bin/g++-9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.o -c /home/xsarm06/Development/cpp_projects/optimization_methods_x64/unit_tests/om_fibonacci_t.cpp

unit_tests/CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.i"
	cd /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests && /bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xsarm06/Development/cpp_projects/optimization_methods_x64/unit_tests/om_fibonacci_t.cpp > CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.i

unit_tests/CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.s"
	cd /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests && /bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xsarm06/Development/cpp_projects/optimization_methods_x64/unit_tests/om_fibonacci_t.cpp -o CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.s

# Object files for target fibonacci_test
fibonacci_test_OBJECTS = \
"CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.o"

# External object files for target fibonacci_test
fibonacci_test_EXTERNAL_OBJECTS =

unit_tests/fibonacci_test: unit_tests/CMakeFiles/fibonacci_test.dir/om_fibonacci_t.cpp.o
unit_tests/fibonacci_test: unit_tests/CMakeFiles/fibonacci_test.dir/build.make
unit_tests/fibonacci_test: liboptimization_methods.a
unit_tests/fibonacci_test: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so.1.71.0
unit_tests/fibonacci_test: unit_tests/CMakeFiles/fibonacci_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fibonacci_test"
	cd /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fibonacci_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unit_tests/CMakeFiles/fibonacci_test.dir/build: unit_tests/fibonacci_test

.PHONY : unit_tests/CMakeFiles/fibonacci_test.dir/build

unit_tests/CMakeFiles/fibonacci_test.dir/clean:
	cd /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests && $(CMAKE_COMMAND) -P CMakeFiles/fibonacci_test.dir/cmake_clean.cmake
.PHONY : unit_tests/CMakeFiles/fibonacci_test.dir/clean

unit_tests/CMakeFiles/fibonacci_test.dir/depend:
	cd /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xsarm06/Development/cpp_projects/optimization_methods_x64 /home/xsarm06/Development/cpp_projects/optimization_methods_x64/unit_tests /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests /home/xsarm06/Development/cpp_projects/optimization_methods_x64/build/unit_tests/CMakeFiles/fibonacci_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unit_tests/CMakeFiles/fibonacci_test.dir/depend

