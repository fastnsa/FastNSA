add_library(fastnsa_capi SHARED fastnsa_capi.cpp)

target_include_directories(fastnsa_capi
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fastnsa>
)

target_compile_definitions(fastnsa_capi PRIVATE FASTNSA_BUILDING_DLL)
target_link_libraries(fastnsa_capi
    PRIVATE
        fastnsa_nne
        fastnsa_stars
        fastnsa_earth
        fastnsa_kdv
)

if(APPLE)
    target_link_options(fastnsa_capi PRIVATE $<$<CONFIG:Release>:-Wl,-dead_strip>)
elseif(MSVC)
    target_link_options(fastnsa_capi PRIVATE $<$<CONFIG:Release>:/OPT:REF> $<$<CONFIG:Release>:/OPT:ICF>)
elseif(MINGW)
    target_link_options(fastnsa_capi PRIVATE $<$<CONFIG:Release>:-Wl,--gc-sections>)
elseif(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_options(fastnsa_capi PRIVATE
        $<$<CONFIG:Release>:-Wl,--gc-sections>
        "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/fastnsa_exports.map"
    )
    set_property(
        TARGET fastnsa_capi
        APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/fastnsa_exports.map"
    )
endif()

set_target_properties(fastnsa_capi PROPERTIES
    OUTPUT_NAME fastnsa
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
)

install(
    TARGETS fastnsa_capi
    LIBRARY DESTINATION "${FASTNSA_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${FASTNSA_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
    FILES fastnsa_capi.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/fastnsa"
)
