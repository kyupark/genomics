# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/el114/kyu/code/mybamtools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/el114/kyu/code/mybamtools/build

# Include any dependencies generated for this target.
include src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/depend.make

# Include the progress variables for this target.
include src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/progress.make

# Include the compile flags for this target's objects.
include src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/flags.make

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/flags.make
src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o: ../src/third_party/jsoncpp/json_reader.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/jsoncppST.dir/json_reader.cpp.o -c /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_reader.cpp

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/jsoncppST.dir/json_reader.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_reader.cpp > CMakeFiles/jsoncppST.dir/json_reader.cpp.i

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/jsoncppST.dir/json_reader.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_reader.cpp -o CMakeFiles/jsoncppST.dir/json_reader.cpp.s

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.requires:
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.requires

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.provides: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.requires
	$(MAKE) -f src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/build.make src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.provides.build
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.provides

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.provides.build: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/flags.make
src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o: ../src/third_party/jsoncpp/json_value.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/jsoncppST.dir/json_value.cpp.o -c /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_value.cpp

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/jsoncppST.dir/json_value.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_value.cpp > CMakeFiles/jsoncppST.dir/json_value.cpp.i

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/jsoncppST.dir/json_value.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_value.cpp -o CMakeFiles/jsoncppST.dir/json_value.cpp.s

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.requires:
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.requires

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.provides: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.requires
	$(MAKE) -f src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/build.make src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.provides.build
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.provides

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.provides.build: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/flags.make
src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o: ../src/third_party/jsoncpp/json_writer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/jsoncppST.dir/json_writer.cpp.o -c /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_writer.cpp

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/jsoncppST.dir/json_writer.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_writer.cpp > CMakeFiles/jsoncppST.dir/json_writer.cpp.i

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/jsoncppST.dir/json_writer.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp/json_writer.cpp -o CMakeFiles/jsoncppST.dir/json_writer.cpp.s

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.requires:
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.requires

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.provides: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.requires
	$(MAKE) -f src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/build.make src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.provides.build
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.provides

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.provides.build: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o

# Object files for target jsoncppST
jsoncppST_OBJECTS = \
"CMakeFiles/jsoncppST.dir/json_reader.cpp.o" \
"CMakeFiles/jsoncppST.dir/json_value.cpp.o" \
"CMakeFiles/jsoncppST.dir/json_writer.cpp.o"

# External object files for target jsoncppST
jsoncppST_EXTERNAL_OBJECTS =

../lib/libjsoncpp.a: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o
../lib/libjsoncpp.a: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o
../lib/libjsoncpp.a: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o
../lib/libjsoncpp.a: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/build.make
../lib/libjsoncpp.a: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../../../lib/libjsoncpp.a"
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && $(CMAKE_COMMAND) -P CMakeFiles/jsoncppST.dir/cmake_clean_target.cmake
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/jsoncppST.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/build: ../lib/libjsoncpp.a
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/build

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/requires: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_reader.cpp.o.requires
src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/requires: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_value.cpp.o.requires
src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/requires: src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/json_writer.cpp.o.requires
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/requires

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/clean:
	cd /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp && $(CMAKE_COMMAND) -P CMakeFiles/jsoncppST.dir/cmake_clean.cmake
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/clean

src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/depend:
	cd /home/el114/kyu/code/mybamtools/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/el114/kyu/code/mybamtools /home/el114/kyu/code/mybamtools/src/third_party/jsoncpp /home/el114/kyu/code/mybamtools/build /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp /home/el114/kyu/code/mybamtools/build/src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/third_party/jsoncpp/CMakeFiles/jsoncppST.dir/depend
