add_executable(integratorTest IntegratorTest.cpp)
add_executable(directsolvertest DirectSolverTest.cpp)

target_link_libraries(integratorTest fastibem hpro)
target_link_libraries(directsolvertest fastibem) 

add_test(integrationTest integratorTest ${fast-ibem_SOURCE_DIR}/input/sphere.dat)
add_test(directsolvertest ${fast-ibem_SOURCE_DIR}/input/sphere.dat)
