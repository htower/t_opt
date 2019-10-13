cmake_minimum_required(VERSION 3.0)
project(t_opt)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -Wpedantic")

# set(CMAKE_C_FLAGS "-O2 -Wall -Wextra")

include_directories(ext)
include_directories(t_opt)

add_executable(t_opt

ext/fmt/src/format.cc
ext/fmt/src/posix.cc

t_opt/core/types.cpp
t_opt/core/problem.cpp
t_opt/core/method.cpp
t_opt/core/logger.cpp

t_opt/utils.cpp

t_opt/line_search/h_simple.cpp
t_opt/line_search/parabolic.cpp

t_opt/local/afgm.cpp
t_opt/local/agmsdr.cpp
t_opt/local/cg.cpp
t_opt/local/fgm.cpp
t_opt/local/gdm.cpp
t_opt/local/lbfgs.cpp
t_opt/local/ufgm.cpp

problems/quadratic.cpp

problems/transport/src/tntp.cpp
problems/transport/src/sdm.cpp
problems/transport/smvsdm2.cpp
problems/transport/tsdm.cpp
problems/transport/lpsdm.cpp

main.cpp)

install(TARGETS t_opt RUNTIME DESTINATION bin)
