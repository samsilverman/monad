cmake_minimum_required(VERSION 3.10)

set(TEST_SOURCES
    grid/hex8_topology_test.cpp
    grid/hex20_topology_test.cpp
    grid/quad4_topology_test.cpp
    grid/quad8_topology_test.cpp
    grid/structured_grid_test.cpp
)

add_executable(monad_tests ${TEST_SOURCES})

target_include_directories(monad_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/tests
)

target_link_libraries(monad_tests PRIVATE
    monad
    Catch2::Catch2WithMain
)

target_compile_options(monad_tests PRIVATE ${MONAD_WARNING_FLAGS})
