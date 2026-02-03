set(MONAD_SOURCES
    detail/mean.cpp
    fem/element/hex8.cpp
    fem/element/hex20.cpp
    fem/element/quad4.cpp
    fem/element/quad8.cpp
    grid/hex8_grid.cpp
    grid/hex20_grid.cpp
    grid/quad4_grid.cpp
    grid/quad8_grid.cpp
    io/gmsh/write_gmsh_header.cpp
    material/mechanical/linear_elastic_2d.cpp
    material/mechanical/linear_elastic_3d.cpp
    solver/solver_options.cpp
)

add_library(monad ${MONAD_SOURCES})

target_include_directories(monad PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(monad PUBLIC cxx_std_17)

target_link_libraries(monad PUBLIC Eigen3::Eigen)

if(MONAD_USE_OPENMP)
    target_link_libraries(monad PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_options(monad PRIVATE ${MONAD_WARNING_FLAGS})
