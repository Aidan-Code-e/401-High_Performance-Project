cmake_minimum_required(VERSION 3.28.1)
project(CREATESIGS LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_FLAGS_DEBUG "-O3 -mavx2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -mavx2")

set(RESOURCE_FILES "small.fasta"
                   "qut2.fasta"
                   "qut3.fasta")

file(COPY ${RESOURCE_FILES}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include(FetchContent)

option(TRACY_ENABLE "" ON)
option(TRACY_NO_EXIT "" ON)

FetchContent_Declare(
     tracy
     GIT_REPOSITORY https://github.com/wolfpld/tracy.git
     GIT_TAG        v0.12.2
     GIT_SHALLOW    TRUE
     GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(tracy)

add_executable(CREATESIGS ISAAC-rand.cpp signature.cpp "${tracy_SOURCE_DIR}/public/TracyClient.cpp")

target_link_libraries(CREATESIGS PUBLIC TracyClient)
target_include_directories(CREATESIGS PUBLIC "${tracy_SOURCE_DIR}/public/tracy/Tracy.hpp")

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
    target_link_libraries(CREATESIGS PUBLIC OpenMP::OpenMP_CXX)
endif()
