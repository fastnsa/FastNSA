
set(NNE_CORE_SOURCES
        core/BoundablePair.cpp
        core/NNE.cpp
        core/SpatialGrid.cpp
        distance/PointLineDistance.cpp
        geometry/distance.cpp
        index/MySTRtree.cpp
)

add_library(fastnsa_nne STATIC
    ${NNE_CORE_SOURCES}
)

target_include_directories(fastnsa_nne
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fastnsa_nne
    PUBLIC
        GEOS::geos
        OpenMP::OpenMP_CXX
)
