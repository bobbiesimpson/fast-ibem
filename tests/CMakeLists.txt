add_executable(integratorTest IntegratorTest.cpp)
add_executable(emagtest EmagTest.cpp)
add_executable(rcstest RCSTest.cpp)
add_executable(projectiontest ProjectionTest.cpp)

target_link_libraries(integratorTest fastibem hpro ${THIS_LINK_LIBS})
target_link_libraries(emagtest fastibem hpro ${THIS_LINK_LIBS})
target_link_libraries(rcstest fastibem hpro ${THIS_LINK_LIBS})
target_link_libraries(projectiontest fastibem hpro ${THIS_LINK_LIBS})

add_test(integrationTest integratorTest ${fast-ibem_SOURCE_DIR}/input/sphere.dat)
add_test(emagtest emagtest)
add_test(rcstest rcstest)
add_test(projectiontest projectiontest)

