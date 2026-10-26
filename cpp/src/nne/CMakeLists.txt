
set(NNE_CORE_SOURCES
        core/NNE.cpp
        core/SpatialGrid.cpp
        core/Threading.cpp
        geometry/Polyline.cpp
        index/PolylineBVH.cpp
)

add_library(fastnsa_nne STATIC
    ${NNE_CORE_SOURCES}
)

target_include_directories(fastnsa_nne
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fastnsa_nne PUBLIC Threads::Threads)
