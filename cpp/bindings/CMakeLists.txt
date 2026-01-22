add_library(fastnsa_capi SHARED
    fastnsa_capi.cpp
)

target_include_directories(fastnsa_capi
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fastnsa_capi
    PRIVATE
        fastnsa_nne
        fastnsa_stars
        OpenMP::OpenMP_CXX
)

set_target_properties(fastnsa_capi PROPERTIES
    OUTPUT_NAME "fastnsa"
)


add_custom_command(TARGET fastnsa_capi POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:fastnsa_capi>
            ${PROJECT_SOURCE_DIR}/../python/fastnsa/bindings/
)
