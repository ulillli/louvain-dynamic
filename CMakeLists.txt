cmake_minimum_required(VERSION 3.15)
project(Louvain CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

file(GLOB SOURCES "src/*.cpp")

add_executable(louvain_app main.cpp ${SOURCES})

target_include_directories(louvain_app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

