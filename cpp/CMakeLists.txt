cmake_minimum_required(VERSION 3.18)

project(fastnsa VERSION 0.3.0 LANGUAGES C CXX)

include(CheckIPOSupported)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Export only the stable C ABI when FastNSA is loaded into a host process.
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(FASTNSA_BUILD_CAPI "Build the FastNSA C API" ON)
option(FASTNSA_BUILD_TESTS "Build FastNSA tests" ON)
option(FASTNSA_ENABLE_IPO "Enable release interprocedural optimization" ON)
set(
    FASTNSA_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}"
    CACHE STRING "Install directory for the FastNSA shared library"
)
# Put functions and data in separate sections so the final shared library can
# discard unused native implementation code.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        $<$<CONFIG:Release>:-ffunction-sections>
        $<$<CONFIG:Release>:-fdata-sections>
    )
elseif(MSVC)
    add_compile_options($<$<CONFIG:Release>:/Gy> $<$<CONFIG:Release>:/Gw>)
endif()

if(FASTNSA_ENABLE_IPO)
    check_ipo_supported(RESULT FASTNSA_IPO_SUPPORTED OUTPUT FASTNSA_IPO_ERROR)
    if(FASTNSA_IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL ON)
    else()
        message(WARNING "IPO/LTO is unavailable: ${FASTNSA_IPO_ERROR}")
    endif()
endif()

find_package(Threads REQUIRED)

add_subdirectory(src/nne)
add_subdirectory(src/stars)
add_subdirectory(src/earth)
add_subdirectory(src/kdv)

if(FASTNSA_BUILD_CAPI)
    add_subdirectory(bindings)
endif()

if(FASTNSA_BUILD_TESTS)
    if(NOT FASTNSA_BUILD_CAPI)
        message(FATAL_ERROR "FASTNSA_BUILD_TESTS requires FASTNSA_BUILD_CAPI=ON")
    endif()
    enable_testing()
    add_subdirectory(tests)
endif()
