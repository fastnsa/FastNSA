function(fastnsa_add_capi_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE fastnsa_capi)
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/bindings)
    if(WIN32)
        add_custom_command(
            TARGET ${name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:fastnsa_capi>
                $<TARGET_FILE_DIR:${name}>
        )
    endif()
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES TIMEOUT 30)
endfunction()

fastnsa_add_capi_test(fastnsa_capi_network test_capi_network.cpp)
fastnsa_add_capi_test(fastnsa_capi_mapping test_capi_mapping.cpp)
fastnsa_add_capi_test(fastnsa_capi_statistics test_capi_statistics.cpp)
fastnsa_add_capi_test(fastnsa_capi_nkdv test_capi_nkdv.cpp)
target_link_libraries(fastnsa_capi_nkdv PRIVATE Threads::Threads)
fastnsa_add_capi_test(fastnsa_capi_kfunction_safety test_capi_kfunction_safety.cpp)
target_link_libraries(fastnsa_capi_kfunction_safety PRIVATE Threads::Threads)
set_tests_properties(fastnsa_capi_kfunction_safety PROPERTIES TIMEOUT 120)
