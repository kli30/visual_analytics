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
include CMakeFiles/triSurface.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/triSurface.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/triSurface.dir/flags.make

CMakeFiles/triSurface.dir/triSurface.cpp.o: CMakeFiles/triSurface.dir/flags.make
CMakeFiles/triSurface.dir/triSurface.cpp.o: ../triSurface.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kaiming/programs/visualAnalytics/kml/src/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/triSurface.dir/triSurface.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/triSurface.dir/triSurface.cpp.o -c /home/kaiming/programs/visualAnalytics/kml/src/triSurface.cpp

CMakeFiles/triSurface.dir/triSurface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/triSurface.dir/triSurface.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/kaiming/programs/visualAnalytics/kml/src/triSurface.cpp > CMakeFiles/triSurface.dir/triSurface.cpp.i

CMakeFiles/triSurface.dir/triSurface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/triSurface.dir/triSurface.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/kaiming/programs/visualAnalytics/kml/src/triSurface.cpp -o CMakeFiles/triSurface.dir/triSurface.cpp.s

CMakeFiles/triSurface.dir/triSurface.cpp.o.requires:
.PHONY : CMakeFiles/triSurface.dir/triSurface.cpp.o.requires

CMakeFiles/triSurface.dir/triSurface.cpp.o.provides: CMakeFiles/triSurface.dir/triSurface.cpp.o.requires
	$(MAKE) -f CMakeFiles/triSurface.dir/build.make CMakeFiles/triSurface.dir/triSurface.cpp.o.provides.build
.PHONY : CMakeFiles/triSurface.dir/triSurface.cpp.o.provides

CMakeFiles/triSurface.dir/triSurface.cpp.o.provides.build: CMakeFiles/triSurface.dir/triSurface.cpp.o

# Object files for target triSurface
triSurface_OBJECTS = \
"CMakeFiles/triSurface.dir/triSurface.cpp.o"

# External object files for target triSurface
triSurface_EXTERNAL_OBJECTS =

libtriSurface.a: CMakeFiles/triSurface.dir/triSurface.cpp.o
libtriSurface.a: CMakeFiles/triSurface.dir/build.make
libtriSurface.a: CMakeFiles/triSurface.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libtriSurface.a"
	$(CMAKE_COMMAND) -P CMakeFiles/triSurface.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/triSurface.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/triSurface.dir/build: libtriSurface.a
.PHONY : CMakeFiles/triSurface.dir/build

CMakeFiles/triSurface.dir/requires: CMakeFiles/triSurface.dir/triSurface.cpp.o.requires
.PHONY : CMakeFiles/triSurface.dir/requires

CMakeFiles/triSurface.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/triSurface.dir/cmake_clean.cmake
.PHONY : CMakeFiles/triSurface.dir/clean

CMakeFiles/triSurface.dir/depend:
	cd /home/kaiming/programs/visualAnalytics/kml/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kaiming/programs/visualAnalytics/kml/src /home/kaiming/programs/visualAnalytics/kml/src /home/kaiming/programs/visualAnalytics/kml/src/build /home/kaiming/programs/visualAnalytics/kml/src/build /home/kaiming/programs/visualAnalytics/kml/src/build/CMakeFiles/triSurface.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/triSurface.dir/depend

