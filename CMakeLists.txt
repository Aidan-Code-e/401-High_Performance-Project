cmake_minimum_required(VERSION 3.10.0)
project(CREATESIGS LANGUAGES C CXX)

set(RESOURCE_FILES "small.fasta"
                   "qut2.fasta"
                   "qut3.fasta")

file(COPY ${RESOURCE_FILES}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(CREATESIGS ISAAC-rand.cpp signature.cpp)
