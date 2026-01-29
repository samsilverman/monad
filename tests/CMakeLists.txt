cmake_minimum_required(VERSION 3.10)

file(GLOB_RECURSE TEST_SOURCES ${CMAKE_SOURCE_DIR}/tests/*.cpp)

add_executable(monad_tests ${TEST_SOURCES})

target_include_directories(monad_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/tests
)

target_link_libraries(monad_tests PRIVATE
    monad
    Catch2::Catch2WithMain
)

target_compile_options(monad_tests PRIVATE ${MONAD_WARNING_FLAGS})
