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
include CMakeFiles/jointModel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/jointModel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/jointModel.dir/flags.make

CMakeFiles/jointModel.dir/jointModel.cpp.o: CMakeFiles/jointModel.dir/flags.make
CMakeFiles/jointModel.dir/jointModel.cpp.o: ../jointModel.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kaiming/programs/visualAnalytics/kml/src/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/jointModel.dir/jointModel.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/jointModel.dir/jointModel.cpp.o -c /home/kaiming/programs/visualAnalytics/kml/src/jointModel.cpp

CMakeFiles/jointModel.dir/jointModel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/jointModel.dir/jointModel.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/kaiming/programs/visualAnalytics/kml/src/jointModel.cpp > CMakeFiles/jointModel.dir/jointModel.cpp.i

CMakeFiles/jointModel.dir/jointModel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/jointModel.dir/jointModel.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/kaiming/programs/visualAnalytics/kml/src/jointModel.cpp -o CMakeFiles/jointModel.dir/jointModel.cpp.s

CMakeFiles/jointModel.dir/jointModel.cpp.o.requires:
.PHONY : CMakeFiles/jointModel.dir/jointModel.cpp.o.requires

CMakeFiles/jointModel.dir/jointModel.cpp.o.provides: CMakeFiles/jointModel.dir/jointModel.cpp.o.requires
	$(MAKE) -f CMakeFiles/jointModel.dir/build.make CMakeFiles/jointModel.dir/jointModel.cpp.o.provides.build
.PHONY : CMakeFiles/jointModel.dir/jointModel.cpp.o.provides

CMakeFiles/jointModel.dir/jointModel.cpp.o.provides.build: CMakeFiles/jointModel.dir/jointModel.cpp.o

# Object files for target jointModel
jointModel_OBJECTS = \
"CMakeFiles/jointModel.dir/jointModel.cpp.o"

# External object files for target jointModel
jointModel_EXTERNAL_OBJECTS =

libjointModel.a: CMakeFiles/jointModel.dir/jointModel.cpp.o
libjointModel.a: CMakeFiles/jointModel.dir/build.make
libjointModel.a: CMakeFiles/jointModel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libjointModel.a"
	$(CMAKE_COMMAND) -P CMakeFiles/jointModel.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/jointModel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/jointModel.dir/build: libjointModel.a
.PHONY : CMakeFiles/jointModel.dir/build

CMakeFiles/jointModel.dir/requires: CMakeFiles/jointModel.dir/jointModel.cpp.o.requires
.PHONY : CMakeFiles/jointModel.dir/requires

CMakeFiles/jointModel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/jointModel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/jointModel.dir/clean

CMakeFiles/jointModel.dir/depend:
	cd /home/kaiming/programs/visualAnalytics/kml/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kaiming/programs/visualAnalytics/kml/src /home/kaiming/programs/visualAnalytics/kml/src /home/kaiming/programs/visualAnalytics/kml/src/build /home/kaiming/programs/visualAnalytics/kml/src/build /home/kaiming/programs/visualAnalytics/kml/src/build/CMakeFiles/jointModel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/jointModel.dir/depend

