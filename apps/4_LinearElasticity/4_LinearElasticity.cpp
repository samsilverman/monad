/**
 * @brief Command-line tool for obtaining the homogenized 2D stiffness tensor for a unit cell.
 *
 * This program obtains the homogenized 2D stiffness tensor for a 15x15 Quad8 grid mesh.
 * The following are saved to Gmsh files:
 *
 * - The density field.
 * - The macroscopic displacements ū₁₁.
 * - The microscopic displacements ũ₁₁.
 * - The total displacements u₁₁=ū₁₁+ũ₁₁.
 *
 * The Gmsh files are written to:
 *
 *      `/path/to/apps/4_LinearElasticity`
 *
 * @param E Young's modulus of the base material (default=1).
 * @param nu Poisson's ratio of the base material (default=0.3).
 *
 * Usage:
 *      $ 4_LinearElasticity [E] [nu]
 *
 * Example:
 *      $ 4_LinearElasticity 1 0.3
 *      ---Homogenized stiffness tensor---
 *          0.131143    0.0567244  8.43549e-15
 *         0.0567244    0.0586181 -3.55764e-15
 *       8.43549e-15 -3.55764e-15    0.0150901
 *      Saved to /path/to/apps/4_LinearElasticity
 */
#include <iostream>
#include <filesystem>
#include <string>
#include "monad/monad.hpp"

int main(int argc, char* argv[]) {
    if (argc < 1 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " [E] [nu]\n";
        return 1;
    }

    double E = 1.0;
    double nu = 0.3;

    try {
        if (argc > 1) {
            E = std::stod(argv[1]);
        }
        if (argc > 2) {
            nu = std::stod(argv[2]);
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "4_LinearElasticity: Invalid argument: " << e.what() << std::endl;
        return 1;
    }

    if (E <= 0) {
        std::cerr << "4_LinearElasticity: E (" << E <<  ") must be positive.\n";
        return 1;
    }
    if (nu <= -1.0 || nu >= 0.5) {
        std::cerr << "4_LinearElasticity: nu (" << nu <<  ") must be in range (-1,0.5).\n";
        return 1;
    }

    monad::Quad8Grid grid({32, 32}, {1.0, 1.0});

    const monad::LinearElasticMaterial2d material(E, nu, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);

    const auto folder = std::filesystem::path(__FILE__).parent_path();
    const auto csvFile = folder / "data.csv";

    grid.setDensitiesFile(csvFile.string());

    monad::SolverOptions options;
    options.fields = monad::FieldSave::All;

    const monad::LinearElasticSolver solver(grid, material);

    const auto results = solver.solve(options);

    std::cout << "---Homogenized stiffness tensor---\n" << results.CBar << std::endl;

    auto file = folder / "density.msh";

    monad::saveGrid(grid, file.string(), true);

    file = folder / "uMacro.msh";

    monad::saveGridAndField(grid, results.uMacro[0], file.string(), "Macro displacement");

    file = folder / "uMicro.msh";

    monad::saveGridAndField(grid, results.uMicro[0], file.string(), "Micro displacement");

    file = folder / "u.msh";

    monad::saveGridAndField(grid, results.u[0], file.string(), "Displacement");

    std::cout << "Saved to " + file.string() << std::endl;

    return 0;
}
