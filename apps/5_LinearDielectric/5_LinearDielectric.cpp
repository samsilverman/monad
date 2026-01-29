/**
 * @brief Command-line tool for obtaining the homogenized 2D permittivity tensor for a unit cell.
 *
 * This program obtains the homogenized 2d permittivity tensor for a 15x15 Quad8 grid mesh.
 * The following are saved to Gmsh files:
 *
 * - The density field.
 * - The macroscopic electric potentials φ̄₁₁
 * - The microscopic electric potentials φ̃₁₁
 * - The total electric potentials φ₁₁=φ̄₁₁+φ̃₁₁
 *
 * The Gmsh files are written to:
 *
 *      `/path/to/apps/5_LinearDielectric`
 *
 * @param epsilon Permittivity constant (default=1).
 *
 * Usage:
 *      $ 5_LinearDielectric [epsilon]
 *
 * Example:
 *      $ 5_LinearDielectric 1
 *      ---Homogenized permittivity tensor---
 *         0.189822 2.5754e-14
 *       2.5754e-14   0.121198
 *      Saved to /path/to/apps/5_LinearDielectric
 */
#include <iostream>
#include <filesystem>
#include <string>
#include "monad/monad.hpp"

int main(int argc, char* argv[]) {
    if (argc < 1 || argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [epsilon]\n";
        return 1;
    }

    double epsilon = 1.0;

    try {
        if (argc > 1) {
            epsilon = std::stod(argv[1]);
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "5_LinearDielectric: Invalid argument: " << e.what() << std::endl;
        return 1;
    }

    if (epsilon <= 0) {
        std::cerr << "5_LinearDielectric: epsilon (" << epsilon <<  ") must be positive.\n";
        return 1;
    }

    monad::Quad8Grid grid({15, 15}, {1.0, 1.0});

    const monad::LinearDielectricMaterial2d material(epsilon);

    const auto folder = std::filesystem::path(__FILE__).parent_path();
    const auto csvFile = folder / "data.csv";

    grid.setDensitiesFile(csvFile.string());

    monad::SolverOptions options;
    options.fields = monad::FieldSave::All;

    const monad::LinearDielectricSolver solver(grid, material);

    // Results follow the generic transport naming:
    // results.KBar → ε̄
    // results.phiMacro → φ̄
    // results.phiMicro → φ̃
    // results.phi → φ
    const auto results = solver.solve(options);

    std::cout << "---Homogenized permittivity tensor---\n" << results.KBar << std::endl;

    auto file = folder / "density.msh";

    monad::saveGrid(grid, file.string(), true);

    file = folder / "phiMacro.msh";

    monad::saveGridAndField(grid, results.phiMacro[0], file.string(), "Macro electric potential");

    file = folder / "phiMicro.msh";

    monad::saveGridAndField(grid, results.phiMicro[0], file.string(), "Micro electric potential");

    file = folder / "phi.msh";

    monad::saveGridAndField(grid, results.phi[0], file.string(), "Electric potential");

    std::cout << "Saved to " + file.string() << std::endl;

    return 0;
}
