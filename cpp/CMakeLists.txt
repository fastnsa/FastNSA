cmake_minimum_required(VERSION 3.10)
project(fastnsa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


option(FASTNSA_BUILD_CAPI "Build C API" ON)


find_package(GEOS REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

add_subdirectory(src/nne)
add_subdirectory(src/stars)


if (FASTNSA_BUILD_CAPI)
    add_subdirectory(bindings)
endif()


