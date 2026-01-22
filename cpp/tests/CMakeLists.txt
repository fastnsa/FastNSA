add_executable(test_fastnsa test_fastnsa.cpp)

target_link_libraries(test_fastnsa 
    PRIVATE 
        fastnsa_capi 
        OpenMP::OpenMP_CXX
)

target_include_directories(test_fastnsa PRIVATE ${CMAKE_SOURCE_DIR}/bindings)


add_test(NAME fastnsa_test COMMAND test_fastnsa)