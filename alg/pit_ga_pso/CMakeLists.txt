cmake_minimum_required(VERSION 3.20)
cmake_policy(VERSION 3.20)

# Set policy to suppress FindBoost deprecation warning
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

# Set vcpkg toolchain file
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../../vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "Vcpkg toolchain file")

# Set vcpkg triplet for MinGW
set(VCPKG_TARGET_TRIPLET "x64-mingw-dynamic" CACHE STRING "")

project(pit_ga_pso LANGUAGES CXX)

# 1) Estándar C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS  OFF)

# 2) nlohmann_json via vcpkg
find_package(nlohmann_json CONFIG REQUIRED)

# 3) Pagmo (try to find, make it optional for now)
find_package(Pagmo QUIET)
if(NOT Pagmo_FOUND)
    message(STATUS "Pagmo not found - building without Pagmo support")
    # Try to find Boost instead for a simpler version
    if(POLICY CMP0167)
        cmake_policy(SET CMP0167 OLD)
    endif()
    find_package(Boost QUIET)
    if(Boost_FOUND)
        message(STATUS "Boost found - can implement basic algorithms manually")
    endif()
endif()

# 4) Ejecutables
if(Pagmo_FOUND)
    add_executable(pit_ga
      src/loader.cpp
      src/fitness.cpp
      src/main_ga.cpp
    )

    add_executable(pit_pso
      src/loader.cpp
      src/fitness.cpp
      src/main_pso.cpp
    )
else()
    # Fallback: simple implementation without Pagmo
    add_executable(pit_simple
      src/loader.cpp
      src/fitness.cpp
      src/main_simple.cpp
    )
endif()

# 5) Includes propios
if(Pagmo_FOUND)
    target_include_directories(pit_ga PRIVATE
      ${CMAKE_SOURCE_DIR}/include
    )

    target_include_directories(pit_pso PRIVATE
      ${CMAKE_SOURCE_DIR}/include
    )
else()
    target_include_directories(pit_simple PRIVATE
      ${CMAKE_SOURCE_DIR}/include
    )
endif()

# 6) Enlazar librerías
if(Pagmo_FOUND)
    target_link_libraries(pit_ga PRIVATE
      nlohmann_json::nlohmann_json
      pagmo
    )

    target_link_libraries(pit_pso PRIVATE
      nlohmann_json::nlohmann_json
      pagmo
    )
else()
    target_link_libraries(pit_simple PRIVATE
      nlohmann_json::nlohmann_json
    )
endif()
