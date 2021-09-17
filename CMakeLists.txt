cmake_minimum_required(VERSION 3.15)
project(bottom_up)

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
add_executable(bottom_up main.cpp Graph_IO.h hcluster.h htree.h DisjointSet.h Utils.h)