cmake_minimum_required(VERSION 3.10)
project(fastnsa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)


option(FASTNSA_BUILD_CAPI "Build C API" ON)
option(FASTNSA_BUILD_TESTS "Build Tests" ON)
# option(FASTNSA_USE_SYSTEM_GEOS "Use system GEOS" OFF)


set(GEOS_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(GEOS_ENABLE_TOOLS OFF CACHE BOOL "" FORCE)
set(GEOS_ENABLE_CAPI OFF CACHE BOOL "" FORCE)
set(GEOS_ENABLE_PYTHON OFF CACHE BOOL "" FORCE)
set(GEOS_BUILD_STATIC ON CACHE BOOL "" FORCE)
set(GEOS_BUILD_SHARED OFF CACHE BOOL "" FORCE)

add_subdirectory(third_party/geos)
set(GEOS_LIB geos)
set(GEOS_INCLUDE ${CMAKE_SOURCE_DIR}/third_party/geos/include)


find_package(OpenMP REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

add_subdirectory(src/nne)
add_subdirectory(src/stars)


if (FASTNSA_BUILD_CAPI)
    add_subdirectory(bindings)
    if (FASTNSA_BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
endif()


