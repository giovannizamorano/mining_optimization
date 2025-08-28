cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.16)
project(pit_ga_pso LANGUAGES CXX)

# 1) Estándar C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS  OFF)

# 2) nlohmann_json via FetchContent
include(FetchContent)
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.2
)
FetchContent_MakeAvailable(nlohmann_json)

# 3) Pagmo
# Si Pagmo tiene CMake config en tu MSYS2, úsalo:
find_package(Pagmo REQUIRED)

# 4) Ejecutable
add_executable(pit_opt
  src/loader.cpp
  src/fitness.cpp
  src/main.cpp
)

# 5) Includes propios
target_include_directories(pit_opt PRIVATE
  ${CMAKE_SOURCE_DIR}/include
)

# 6) Enlazar librerías
target_link_libraries(pit_opt PRIVATE
  pagmo
  nlohmann_json::nlohmann_json
)
