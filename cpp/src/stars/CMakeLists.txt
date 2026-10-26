
set(STARS_CORE_SOURCES
    count.cpp
    K_function.cpp
    Network.cpp
    sam_count.cpp
    SP.cpp
)

add_library(fastnsa_stars STATIC
    ${STARS_CORE_SOURCES}
)

target_include_directories(fastnsa_stars
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)
