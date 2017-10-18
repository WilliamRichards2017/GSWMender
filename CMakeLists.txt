cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

Project(GSWMender)
set(CMAKE_MACOSX_RPATH 1)
set (CMAKE_CXX_STANDARD 11)

add_library(ArrayUtil SHARED ArrayUtil.cpp)
target_link_libraries(ArrayUtil)

add_library(Class-Node SHARED Class-Node.cpp)
target_link_libraries(Class-Node)

add_library(Class-Alignment SHARED Class-Alignment.cpp)
target_link_libraries(Class-Alignment)

add_library(Class-GraphAlignment SHARED Class-GraphAlignment.cpp)
target_link_libraries(Class-GraphAlignment Class-Alignment Class-Node)

add_library(Traceback SHARED Traceback.cpp)
target_link_libraries(Traceback Class-Node Class-GraphAlignment ArrayUtil)

add_library(PileUp SHARED PileUp.cpp)
target_link_libraries(PileUp Traceback ArrayUtil)

add_executable(gsw gsw.cpp)
target_link_libraries(gsw Class-Node Class-Alignment Class-GraphAlignment Traceback PileUp ArrayUtil)

add_executable(tests tests.cpp)
target_link_libraries(tests Class-Node Class-Alignment Class-GraphAlignment Traceback PileUp)
