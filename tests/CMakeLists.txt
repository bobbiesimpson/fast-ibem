add_executable(integratorTest IntegratorTest.cpp)
add_executable(emagtest EmagTest.cpp)

target_link_libraries(integratorTest fastibem hpro ${THIS_LINK_LIBS})
target_link_libraries(emagtest fastibem hpro ${THIS_LINK_LIBS})

add_test(integrationTest integratorTest ${fast-ibem_SOURCE_DIR}/input/sphere.dat)
add_test(emagtest emagtest)
