
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

add_dependencies(fastnsa_nne geos)


target_include_directories(fastnsa_nne
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${GEOS_INCLUDE}
)

target_link_libraries(fastnsa_nne
    PUBLIC
        ${CMAKE_BINARY_DIR}/third_party/geos/lib/libgeos.a
        OpenMP::OpenMP_CXX
)
