cmake_minimum_required(VERSION 3.10)

project(Matrix_Exp)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${PROJECT_SOURCE_DIR})
add_executable(main main.cpp Lei.cpp  original.cpp  Our.cpp  pre_op.cpp)
