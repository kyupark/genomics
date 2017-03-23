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
include src/utils/CMakeFiles/BamTools-utils.dir/depend.make

# Include the progress variables for this target.
include src/utils/CMakeFiles/BamTools-utils.dir/progress.make

# Include the compile flags for this target's objects.
include src/utils/CMakeFiles/BamTools-utils.dir/flags.make

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o: src/utils/CMakeFiles/BamTools-utils.dir/flags.make
src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o: ../src/utils/bamtools_fasta.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o -c /home/el114/kyu/code/mybamtools/src/utils/bamtools_fasta.cpp

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/utils/bamtools_fasta.cpp > CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.i

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/utils/bamtools_fasta.cpp -o CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.s

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.requires:
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.requires

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.provides: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.requires
	$(MAKE) -f src/utils/CMakeFiles/BamTools-utils.dir/build.make src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.provides.build
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.provides

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.provides.build: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o: src/utils/CMakeFiles/BamTools-utils.dir/flags.make
src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o: ../src/utils/bamtools_options.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o -c /home/el114/kyu/code/mybamtools/src/utils/bamtools_options.cpp

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/utils/bamtools_options.cpp > CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.i

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/utils/bamtools_options.cpp -o CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.s

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.requires:
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.requires

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.provides: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.requires
	$(MAKE) -f src/utils/CMakeFiles/BamTools-utils.dir/build.make src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.provides.build
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.provides

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.provides.build: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o: src/utils/CMakeFiles/BamTools-utils.dir/flags.make
src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o: ../src/utils/bamtools_pileup_engine.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o -c /home/el114/kyu/code/mybamtools/src/utils/bamtools_pileup_engine.cpp

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/utils/bamtools_pileup_engine.cpp > CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.i

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/utils/bamtools_pileup_engine.cpp -o CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.s

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.requires:
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.requires

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.provides: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.requires
	$(MAKE) -f src/utils/CMakeFiles/BamTools-utils.dir/build.make src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.provides.build
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.provides

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.provides.build: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o: src/utils/CMakeFiles/BamTools-utils.dir/flags.make
src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o: ../src/utils/bamtools_utilities.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o -c /home/el114/kyu/code/mybamtools/src/utils/bamtools_utilities.cpp

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/utils/bamtools_utilities.cpp > CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.i

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/utils/bamtools_utilities.cpp -o CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.s

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.requires:
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.requires

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.provides: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.requires
	$(MAKE) -f src/utils/CMakeFiles/BamTools-utils.dir/build.make src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.provides.build
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.provides

src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.provides.build: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o

# Object files for target BamTools-utils
BamTools__utils_OBJECTS = \
"CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o" \
"CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o" \
"CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o" \
"CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o"

# External object files for target BamTools-utils
BamTools__utils_EXTERNAL_OBJECTS =

../lib/libbamtools-utils.so.0.9.0: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o
../lib/libbamtools-utils.so.0.9.0: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o
../lib/libbamtools-utils.so.0.9.0: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o
../lib/libbamtools-utils.so.0.9.0: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o
../lib/libbamtools-utils.so.0.9.0: src/utils/CMakeFiles/BamTools-utils.dir/build.make
../lib/libbamtools-utils.so.0.9.0: ../lib/libbamtools.so.0.9.0
../lib/libbamtools-utils.so.0.9.0: src/utils/CMakeFiles/BamTools-utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../../../lib/libbamtools-utils.so"
	cd /home/el114/kyu/code/mybamtools/build/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BamTools-utils.dir/link.txt --verbose=$(VERBOSE)
	cd /home/el114/kyu/code/mybamtools/build/src/utils && $(CMAKE_COMMAND) -E cmake_symlink_library ../../../lib/libbamtools-utils.so.0.9.0 ../../../lib/libbamtools-utils.so.0.9.0 ../../../lib/libbamtools-utils.so

../lib/libbamtools-utils.so: ../lib/libbamtools-utils.so.0.9.0

# Rule to build all files generated by this target.
src/utils/CMakeFiles/BamTools-utils.dir/build: ../lib/libbamtools-utils.so
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/build

src/utils/CMakeFiles/BamTools-utils.dir/requires: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_fasta.cpp.o.requires
src/utils/CMakeFiles/BamTools-utils.dir/requires: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_options.cpp.o.requires
src/utils/CMakeFiles/BamTools-utils.dir/requires: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_pileup_engine.cpp.o.requires
src/utils/CMakeFiles/BamTools-utils.dir/requires: src/utils/CMakeFiles/BamTools-utils.dir/bamtools_utilities.cpp.o.requires
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/requires

src/utils/CMakeFiles/BamTools-utils.dir/clean:
	cd /home/el114/kyu/code/mybamtools/build/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/BamTools-utils.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/clean

src/utils/CMakeFiles/BamTools-utils.dir/depend:
	cd /home/el114/kyu/code/mybamtools/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/el114/kyu/code/mybamtools /home/el114/kyu/code/mybamtools/src/utils /home/el114/kyu/code/mybamtools/build /home/el114/kyu/code/mybamtools/build/src/utils /home/el114/kyu/code/mybamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/BamTools-utils.dir/depend

