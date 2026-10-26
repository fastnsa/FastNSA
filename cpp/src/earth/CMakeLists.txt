
set(EARTH_CORE_SOURCES
    count.cpp
    K_function.cpp
    Network.cpp
    range_tree.cpp
    SP.cpp
)

add_library(fastnsa_earth STATIC
    ${EARTH_CORE_SOURCES}
)

target_include_directories(fastnsa_earth
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)
