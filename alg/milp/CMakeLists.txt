cmake_minimum_required(VERSION 3.16)
project(upl_milp_single_period LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -Wpedantic)

# Requiere GUROBI_HOME configurado
if(NOT DEFINED ENV{GUROBI_HOME})
  message(FATAL_ERROR "GUROBI_HOME no definido. Ej.: C:/gurobi1100/win64 o /opt/gurobi1100/linux64")
endif()
set(GUROBI_HOME $ENV{GUROBI_HOME})

include_directories(${GUROBI_HOME}/include)
link_directories(${GUROBI_HOME}/lib)

add_executable(upl_milp_single_period src/main.cpp)

# Nota: en Linux suele ser 'gurobi_c++' y 'gurobi110' (o 'gurobi100', etc.).
# En Windows, el nombre varía por toolset (p.ej., gurobi_c++md2019). Si falla,
# ajusta el nombre exacto según los .lib en ${GUROBI_HOME}/lib.
target_link_libraries(upl_milp_single_period gurobi_c++ gurobi)
