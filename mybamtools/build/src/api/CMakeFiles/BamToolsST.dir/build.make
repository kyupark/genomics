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
include src/api/CMakeFiles/BamToolsST.dir/depend.make

# Include the progress variables for this target.
include src/api/CMakeFiles/BamToolsST.dir/progress.make

# Include the compile flags for this target's objects.
include src/api/CMakeFiles/BamToolsST.dir/flags.make

src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o: ../src/api/BamAlignment.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/BamAlignment.cpp

src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/BamAlignment.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/BamAlignment.cpp > CMakeFiles/BamToolsST.dir/BamAlignment.cpp.i

src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/BamAlignment.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/BamAlignment.cpp -o CMakeFiles/BamToolsST.dir/BamAlignment.cpp.s

src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o

src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o: ../src/api/BamIndex.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/BamIndex.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/BamIndex.cpp

src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/BamIndex.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/BamIndex.cpp > CMakeFiles/BamToolsST.dir/BamIndex.cpp.i

src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/BamIndex.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/BamIndex.cpp -o CMakeFiles/BamToolsST.dir/BamIndex.cpp.s

src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o

src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o: ../src/api/BamMultiReader.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/BamMultiReader.cpp

src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/BamMultiReader.cpp > CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.i

src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/BamMultiReader.cpp -o CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.s

src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o

src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o: ../src/api/BamReader.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/BamReader.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/BamReader.cpp

src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/BamReader.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/BamReader.cpp > CMakeFiles/BamToolsST.dir/BamReader.cpp.i

src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/BamReader.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/BamReader.cpp -o CMakeFiles/BamToolsST.dir/BamReader.cpp.s

src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o

src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o: ../src/api/BamWriter.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/BamWriter.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/BamWriter.cpp

src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/BamWriter.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/BamWriter.cpp > CMakeFiles/BamToolsST.dir/BamWriter.cpp.i

src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/BamWriter.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/BamWriter.cpp -o CMakeFiles/BamToolsST.dir/BamWriter.cpp.s

src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o

src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o: ../src/api/BGZF.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/BGZF.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/BGZF.cpp

src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/BGZF.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/BGZF.cpp > CMakeFiles/BamToolsST.dir/BGZF.cpp.i

src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/BGZF.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/BGZF.cpp -o CMakeFiles/BamToolsST.dir/BGZF.cpp.s

src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o

src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o: ../src/api/internal/BamReader_p.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/internal/BamReader_p.cpp

src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/internal/BamReader_p.cpp > CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.i

src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/internal/BamReader_p.cpp -o CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.s

src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o

src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o: ../src/api/internal/BamStandardIndex_p.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/internal/BamStandardIndex_p.cpp

src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/internal/BamStandardIndex_p.cpp > CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.i

src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/internal/BamStandardIndex_p.cpp -o CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.s

src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o

src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o: ../src/api/internal/BamToolsIndex_p.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/internal/BamToolsIndex_p.cpp

src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/internal/BamToolsIndex_p.cpp > CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.i

src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/internal/BamToolsIndex_p.cpp -o CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.s

src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o

src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o: src/api/CMakeFiles/BamToolsST.dir/flags.make
src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o: ../src/api/internal/BamWriter_p.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/el114/kyu/code/mybamtools/build/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o -c /home/el114/kyu/code/mybamtools/src/api/internal/BamWriter_p.cpp

src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.i"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/el114/kyu/code/mybamtools/src/api/internal/BamWriter_p.cpp > CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.i

src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.s"
	cd /home/el114/kyu/code/mybamtools/build/src/api && /opt/gcc/4.8.5/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/el114/kyu/code/mybamtools/src/api/internal/BamWriter_p.cpp -o CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.s

src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.requires:
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.requires

src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.provides: src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.requires
	$(MAKE) -f src/api/CMakeFiles/BamToolsST.dir/build.make src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.provides.build
.PHONY : src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.provides

src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.provides.build: src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o

# Object files for target BamToolsST
BamToolsST_OBJECTS = \
"CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o" \
"CMakeFiles/BamToolsST.dir/BamIndex.cpp.o" \
"CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o" \
"CMakeFiles/BamToolsST.dir/BamReader.cpp.o" \
"CMakeFiles/BamToolsST.dir/BamWriter.cpp.o" \
"CMakeFiles/BamToolsST.dir/BGZF.cpp.o" \
"CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o" \
"CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o" \
"CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o" \
"CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o"

# External object files for target BamToolsST
BamToolsST_EXTERNAL_OBJECTS =

../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/build.make
../lib/libbamtools.a: src/api/CMakeFiles/BamToolsST.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../../../lib/libbamtools.a"
	cd /home/el114/kyu/code/mybamtools/build/src/api && $(CMAKE_COMMAND) -P CMakeFiles/BamToolsST.dir/cmake_clean_target.cmake
	cd /home/el114/kyu/code/mybamtools/build/src/api && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BamToolsST.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/api/CMakeFiles/BamToolsST.dir/build: ../lib/libbamtools.a
.PHONY : src/api/CMakeFiles/BamToolsST.dir/build

src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/BamAlignment.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/BamIndex.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/BamMultiReader.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/BamReader.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/BamWriter.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/BGZF.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/internal/BamReader_p.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/internal/BamStandardIndex_p.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/internal/BamToolsIndex_p.cpp.o.requires
src/api/CMakeFiles/BamToolsST.dir/requires: src/api/CMakeFiles/BamToolsST.dir/internal/BamWriter_p.cpp.o.requires
.PHONY : src/api/CMakeFiles/BamToolsST.dir/requires

src/api/CMakeFiles/BamToolsST.dir/clean:
	cd /home/el114/kyu/code/mybamtools/build/src/api && $(CMAKE_COMMAND) -P CMakeFiles/BamToolsST.dir/cmake_clean.cmake
.PHONY : src/api/CMakeFiles/BamToolsST.dir/clean

src/api/CMakeFiles/BamToolsST.dir/depend:
	cd /home/el114/kyu/code/mybamtools/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/el114/kyu/code/mybamtools /home/el114/kyu/code/mybamtools/src/api /home/el114/kyu/code/mybamtools/build /home/el114/kyu/code/mybamtools/build/src/api /home/el114/kyu/code/mybamtools/build/src/api/CMakeFiles/BamToolsST.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/api/CMakeFiles/BamToolsST.dir/depend

