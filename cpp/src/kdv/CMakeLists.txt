add_library(fastnsa_kdv STATIC
    alg_NKDV.cpp
    KAF.cpp
    lion.cpp
    shortest_path.cpp
    compute.cpp
)

target_include_directories(fastnsa_kdv
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)
