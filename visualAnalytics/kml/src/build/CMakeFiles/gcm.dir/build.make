# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kaiming/programs/visualAnalytics/kml/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kaiming/programs/visualAnalytics/kml/src/build

# Include any dependencies generated for this target.
include CMakeFiles/gcm.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/gcm.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gcm.dir/flags.make

CMakeFiles/gcm.dir/grangerCausality.cpp.o: CMakeFiles/gcm.dir/flags.make
CMakeFiles/gcm.dir/grangerCausality.cpp.o: ../grangerCausality.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kaiming/programs/visualAnalytics/kml/src/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/gcm.dir/grangerCausality.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/gcm.dir/grangerCausality.cpp.o -c /home/kaiming/programs/visualAnalytics/kml/src/grangerCausality.cpp

CMakeFiles/gcm.dir/grangerCausality.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gcm.dir/grangerCausality.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/kaiming/programs/visualAnalytics/kml/src/grangerCausality.cpp > CMakeFiles/gcm.dir/grangerCausality.cpp.i

CMakeFiles/gcm.dir/grangerCausality.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gcm.dir/grangerCausality.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/kaiming/programs/visualAnalytics/kml/src/grangerCausality.cpp -o CMakeFiles/gcm.dir/grangerCausality.cpp.s

CMakeFiles/gcm.dir/grangerCausality.cpp.o.requires:
.PHONY : CMakeFiles/gcm.dir/grangerCausality.cpp.o.requires

CMakeFiles/gcm.dir/grangerCausality.cpp.o.provides: CMakeFiles/gcm.dir/grangerCausality.cpp.o.requires
	$(MAKE) -f CMakeFiles/gcm.dir/build.make CMakeFiles/gcm.dir/grangerCausality.cpp.o.provides.build
.PHONY : CMakeFiles/gcm.dir/grangerCausality.cpp.o.provides

CMakeFiles/gcm.dir/grangerCausality.cpp.o.provides.build: CMakeFiles/gcm.dir/grangerCausality.cpp.o

# Object files for target gcm
gcm_OBJECTS = \
"CMakeFiles/gcm.dir/grangerCausality.cpp.o"

# External object files for target gcm
gcm_EXTERNAL_OBJECTS =

libgcm.a: CMakeFiles/gcm.dir/grangerCausality.cpp.o
libgcm.a: CMakeFiles/gcm.dir/build.make
libgcm.a: CMakeFiles/gcm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libgcm.a"
	$(CMAKE_COMMAND) -P CMakeFiles/gcm.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gcm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gcm.dir/build: libgcm.a
.PHONY : CMakeFiles/gcm.dir/build

CMakeFiles/gcm.dir/requires: CMakeFiles/gcm.dir/grangerCausality.cpp.o.requires
.PHONY : CMakeFiles/gcm.dir/requires

CMakeFiles/gcm.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gcm.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gcm.dir/clean

CMakeFiles/gcm.dir/depend:
	cd /home/kaiming/programs/visualAnalytics/kml/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kaiming/programs/visualAnalytics/kml/src /home/kaiming/programs/visualAnalytics/kml/src /home/kaiming/programs/visualAnalytics/kml/src/build /home/kaiming/programs/visualAnalytics/kml/src/build /home/kaiming/programs/visualAnalytics/kml/src/build/CMakeFiles/gcm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gcm.dir/depend
