cmake_minimum_required(VERSION 3.15)
project(milp)

set(CMAKE_CXX_STANDARD 17)

# 1) Indicar el toolchain de vcpkg:
set(CMAKE_TOOLCHAIN_FILE
    "C:/Users/Giova/OneDrive/Escritorio/Mining/vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "Vcpkg toolchain file")

# 2) Encuentra la carpeta de include de mlpack (header-only):
find_path(MLPACK_INCLUDE_DIRS 
  NAMES mlpack/core.hpp 
  HINTS "${_VCPKG_INSTALLED_DIR}/include"
)

if (NOT MLPACK_INCLUDE_DIRS)
  message(FATAL_ERROR "No pude encontrar mlpack.hpp. ¿Está instalado?")
endif()

message(STATUS "mlpack include directory: ${MLPACK_INCLUDE_DIRS}")

# 3) Tu ejecutable
add_executable(milp milp_rf.cpp)

# 4) Indica a CMake dónde buscar los headers de mlpack
target_include_directories(milp PRIVATE ${MLPACK_INCLUDE_DIRS})
