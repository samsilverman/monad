/**
 * @brief Command-line tool for obtaining the homogenized 2D stiffness, permittivity, and piezoelectric tensors for a unit cell.
 *
 * This program obtains the homogenized 2d stiffness, permittivity, and piezoelectric tensors for a 15x15 Quad8 grid mesh.
 * The following are saved to Gmsh files:
 *
 * - The density field.
 * - The macroscopic displacements ū₁₁.
 * - The microscopic displacements ũ₁₁.
 * - The total displacements u₁₁=ū₁₁+ũ₁₁.
 * - The macroscopic electric potentials φ̄₁₁
 * - The microscopic electric potentials φ̃₁₁
 * - The total electric potentials φ₁₁=φ̄₁₁+φ̃₁₁
 *
 * The base material's piezoelectric coupling tensor is set to
 *
 * ⎡d  0  0⎤
 * ⎣0 -d -d⎦
 *
 * where d=E/10.
 *
 * The Gmsh files are written to:
 *
 *      `/path/to/apps/6_LinearPiezoelectricity`
 *
 * @param E Young's modulus of the base material (default=1).
 * @param nu Poisson's ratio of the base material (default=0.3).
 * @param epsilon Permittivity constant (default=1).
 *
 * Usage:
 *      $ 6_LinearPiezoelectricity [E] [nu] [epsilon]
 *
 * Example:
 *      $ 6_LinearPiezoelectricity 1 0.3 1
 *      ---Homogenized stiffness tensor---
 *         0.131468   0.0569734 3.09494e-05
 *        0.0569734   0.0588763 3.35592e-05
 *      3.09494e-05 3.35592e-05   0.0151413
 *
 *      ---Homogenized permittivity tensor---
 *          0.190909 -0.000693838
 *      -0.000693838     0.123736
 *
 *      ---Homogenized piezoelectric tensor---
 *        0.0156172  0.00799243 0.000359505
 *       0.00231767   0.0028842  0.00253854
 *      Saved to /path/to/apps/6_LinearPiezoelectricity
 */
#include <iostream>
#include <filesystem>
#include <string>
#include <Eigen/Core>
#include "monad/monad.hpp"

int main(int argc, char* argv[]) {
    if (argc < 1 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " [E] [nu] [epsilon]\n";
        return 1;
    }

    double E = 1.0;
    double nu = 0.3;
    double epsilon = 1.0;

    try {
        if (argc > 1) {
            E = std::stod(argv[1]);
        }
        if (argc > 2) {
            nu = std::stod(argv[2]);
        }
        if (argc > 3) {
            epsilon = std::stod(argv[3]);
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "6_LinearPiezoelectricity: Invalid argument: " << e.what() << std::endl;
        return 1;
    }

    if (E <= 0) {
        std::cerr << "6_LinearPiezoelectricity: E (" << E <<  ") must be positive.\n";
        return 1;
    }
    if (nu <= -1.0 || nu >= 0.5) {
        std::cerr << "6_LinearPiezoelectricity: nu (" << nu <<  ") must be in range (-1,0.5).\n";
        return 1;
    }
    if (epsilon <= 0) {
        std::cerr << "6_LinearPiezoelectricity: epsilon (" << epsilon <<  ") must be positive.\n";
        return 1;
    }

    monad::Quad8Grid grid({15, 15}, {1.0, 1.0});

    const monad::LinearElasticMaterial2d elasticMaterial(E, nu, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const monad::LinearDielectricMaterial2d dielectricMaterial(epsilon);

    const Eigen::Matrix<double, 2, 3> d {
        {E / 10.0, 0.0, 0.0},
        {0.0, E / 10.0, E / 10.0}
    };

    const monad::LinearPiezoelectricMaterial2d material(elasticMaterial, dielectricMaterial, d);

    const auto folder = std::filesystem::path(__FILE__).parent_path();
    const auto csvFile = folder / "data.csv";

    grid.setDensitiesFile(csvFile.string());

    monad::SolverOptions options;
    options.fields = monad::FieldSave::All;

    const monad::LinearPiezoelectricSolver solver(grid, material);
    const auto results = solver.solve(options);

    std::cout << "---Homogenized stiffness tensor---\n" << results.cBar << std::endl;
    std::cout << "\n---Homogenized permittivity tensor---\n" << results.epsilonBar << std::endl;
    std::cout << "\n---Homogenized piezoelectric tensor---\n" << results.dBar << std::endl;

    auto file = folder / "density.msh";

    monad::saveGrid(grid, file.string(), true);

    file = folder / "uMacro.msh";

    monad::saveGridAndField(grid, results.uMacro[0], file.string(), "Macro displacement");

    file = folder / "uMicro.msh";

    monad::saveGridAndField(grid, results.uMicro[0], file.string(), "Micro displacement");

    file = folder / "u.msh";

    monad::saveGridAndField(grid, results.u[0], file.string(), "Displacement");

    file = folder / "phiMacro.msh";

    monad::saveGridAndField(grid, results.phiMacro[0], file.string(), "Macro electric potential");

    file = folder / "phiMicro.msh";

    monad::saveGridAndField(grid, results.phiMicro[0], file.string(), "Micro electric potential");

    file = folder / "phi.msh";

    monad::saveGridAndField(grid, results.phi[0], file.string(), "Electric potential");

    std::cout << "Saved to " + file.string() << std::endl;

    return 0;
}
