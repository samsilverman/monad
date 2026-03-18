cmake_minimum_required(VERSION 3.10)

set(TEST_SOURCES
    detail/eigen_utils_test.cpp
    detail/mean_test.cpp
    fem/element/element_test.cpp
    fem/kernel/mechanical/linear_elastic_kernel_test.cpp
    fem/kernel/multiphysics/linear_piezoelectric_kernel_test.cpp
    fem/kernel/scalar/linear_scalar_diffusion_kernel_test.cpp
    fem/operator/mechanical/linear_elastic_dof_traits_test.cpp
    fem/operator/multiphysics/linear_piezoelectric_dof_traits_test.cpp
    fem/operator/scalar/linear_scalar_diffusion_dof_traits_test.cpp
    field/density_field_test.cpp
    field/make_density_field_from_csv_test.cpp
    field/make_density_field_from_function_test.cpp
    grid/hex8_topology_test.cpp
    grid/hex20_topology_test.cpp
    grid/quad4_topology_test.cpp
    grid/quad8_topology_test.cpp
    grid/structured_grid_test.cpp
    integration/integrate_matrix_test.cpp
    integration/integrate_scalar_test.cpp
    integration/quadrature_rule_test.cpp
    io/gmsh/write_gmsh_densities_test.cpp
    io/gmsh/write_gmsh_elements_test.cpp
    io/gmsh/write_gmsh_header_test.cpp
    io/gmsh/write_gmsh_nodal_field_test.cpp
    io/gmsh/write_gmsh_nodes_test.cpp
    io/save_grid_and_density_field_test.cpp
    io/save_grid_and_nodal_field_test.cpp
    io/save_grid_test.cpp
    material/mechanical/linear_elastic_material_2d_test.cpp
    material/mechanical/linear_elastic_material_3d_test.cpp
    material/multiphysics/linear_piezoelectric_material_2d_test.cpp
    material/multiphysics/linear_piezoelectric_material_3d_test.cpp
    material/transport/linear_transport_material_test.cpp
    material/bounds_test.cpp
    solver/mechanical/linear_elastic_solver_test.cpp
    solver/multiphysics/linear_piezoelectric_solver_test.cpp
    solver/scalar/linear_scalar_diffusion_solver_test.cpp
    solver/solver_options_test.cpp
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
