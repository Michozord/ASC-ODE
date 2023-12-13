cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# use linalg from NGSolve, otherwise comment it out
# find_package(NGSolve CONFIG REQUIRED
#   HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
# )


include_directories(src)
add_subdirectory (src)

# add_executable (test_ode demos/test_ode.cc)

# add_executable (test_newmark demos/test_newmark.cc)

# add_executable (test_alpha demos/test_alpha.cc)

add_executable (test_RK demos/test_RK.cc)

#add_subdirectory (mass_spring)
