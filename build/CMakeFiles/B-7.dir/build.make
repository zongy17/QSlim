# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = D:\cmake-3.21.2-windows-x86_64\bin\cmake.exe

# The command to remove a file.
RM = D:\cmake-3.21.2-windows-x86_64\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\2021Fall\ComputerGraphics\homework\B-7

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\2021Fall\ComputerGraphics\homework\B-7\build

# Include any dependencies generated for this target.
include CMakeFiles/B-7.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/B-7.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/B-7.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/B-7.dir/flags.make

CMakeFiles/B-7.dir/src/glad.c.obj: CMakeFiles/B-7.dir/flags.make
CMakeFiles/B-7.dir/src/glad.c.obj: CMakeFiles/B-7.dir/includes_C.rsp
CMakeFiles/B-7.dir/src/glad.c.obj: ../src/glad.c
CMakeFiles/B-7.dir/src/glad.c.obj: CMakeFiles/B-7.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\2021Fall\ComputerGraphics\homework\B-7\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/B-7.dir/src/glad.c.obj"
	D:\mingw64\bin\x86_64-w64-mingw32-gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/B-7.dir/src/glad.c.obj -MF CMakeFiles\B-7.dir\src\glad.c.obj.d -o CMakeFiles\B-7.dir\src\glad.c.obj -c D:\2021Fall\ComputerGraphics\homework\B-7\src\glad.c

CMakeFiles/B-7.dir/src/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/B-7.dir/src/glad.c.i"
	D:\mingw64\bin\x86_64-w64-mingw32-gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E D:\2021Fall\ComputerGraphics\homework\B-7\src\glad.c > CMakeFiles\B-7.dir\src\glad.c.i

CMakeFiles/B-7.dir/src/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/B-7.dir/src/glad.c.s"
	D:\mingw64\bin\x86_64-w64-mingw32-gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S D:\2021Fall\ComputerGraphics\homework\B-7\src\glad.c -o CMakeFiles\B-7.dir\src\glad.c.s

CMakeFiles/B-7.dir/src/main.cpp.obj: CMakeFiles/B-7.dir/flags.make
CMakeFiles/B-7.dir/src/main.cpp.obj: CMakeFiles/B-7.dir/includes_CXX.rsp
CMakeFiles/B-7.dir/src/main.cpp.obj: ../src/main.cpp
CMakeFiles/B-7.dir/src/main.cpp.obj: CMakeFiles/B-7.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\2021Fall\ComputerGraphics\homework\B-7\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/B-7.dir/src/main.cpp.obj"
	D:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/B-7.dir/src/main.cpp.obj -MF CMakeFiles\B-7.dir\src\main.cpp.obj.d -o CMakeFiles\B-7.dir\src\main.cpp.obj -c D:\2021Fall\ComputerGraphics\homework\B-7\src\main.cpp

CMakeFiles/B-7.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/B-7.dir/src/main.cpp.i"
	D:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\2021Fall\ComputerGraphics\homework\B-7\src\main.cpp > CMakeFiles\B-7.dir\src\main.cpp.i

CMakeFiles/B-7.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/B-7.dir/src/main.cpp.s"
	D:\mingw64\bin\x86_64-w64-mingw32-g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\2021Fall\ComputerGraphics\homework\B-7\src\main.cpp -o CMakeFiles\B-7.dir\src\main.cpp.s

# Object files for target B-7
B__7_OBJECTS = \
"CMakeFiles/B-7.dir/src/glad.c.obj" \
"CMakeFiles/B-7.dir/src/main.cpp.obj"

# External object files for target B-7
B__7_EXTERNAL_OBJECTS =

B-7.exe: CMakeFiles/B-7.dir/src/glad.c.obj
B-7.exe: CMakeFiles/B-7.dir/src/main.cpp.obj
B-7.exe: CMakeFiles/B-7.dir/build.make
B-7.exe: CMakeFiles/B-7.dir/linklibs.rsp
B-7.exe: CMakeFiles/B-7.dir/objects1.rsp
B-7.exe: CMakeFiles/B-7.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\2021Fall\ComputerGraphics\homework\B-7\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable B-7.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\B-7.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/B-7.dir/build: B-7.exe
.PHONY : CMakeFiles/B-7.dir/build

CMakeFiles/B-7.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\B-7.dir\cmake_clean.cmake
.PHONY : CMakeFiles/B-7.dir/clean

CMakeFiles/B-7.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\2021Fall\ComputerGraphics\homework\B-7 D:\2021Fall\ComputerGraphics\homework\B-7 D:\2021Fall\ComputerGraphics\homework\B-7\build D:\2021Fall\ComputerGraphics\homework\B-7\build D:\2021Fall\ComputerGraphics\homework\B-7\build\CMakeFiles\B-7.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/B-7.dir/depend

