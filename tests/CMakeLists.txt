add_executable(integratorTest IntegratorTest.cpp)

target_link_libraries(integratorTest fastibem)

add_test(integrationTest integratorTest ${fast-ibem_SOURCE_DIR}/input/sphere.dat)
