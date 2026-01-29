/**
 * @brief Command-line tool for generating a density field from a function.
 *
 * This program fills the density field of a 32x32x32 Hex8 grid mesh with
 * densities obtained from the continuous function
 *
 * f(x,y,z)=⅓(sin²(10x)+sin²(10y)+sin²(10z))
 *
 * and saves it to a Gmsh-compatible file.
 *
 * The Gmsh file is written to:
 *
 *      `/path/to/apps/3_DensityFunction/output.gmsh`
 *
 * Usage:
 *      $ 3_DensityFunction
 *
 * Example:
 *      $ 3_DensityFunction
 *      Saved to /path/to/apps/3_DensityFunction/output.gmsh
 */
#include <iostream>
#include <filesystem>
#include <cmath>
#include <Eigen/Core>
#include "monad/monad.hpp"

int main() {
    monad::Hex8Grid grid({32, 32, 32}, {1.0, 1.0, 1.0});

    // f(x,y,z)=(sin²(10x)+sin²(10y)+sin²(10z))/3
    auto f = [](const Eigen::Vector3d &point) -> double {
        const double x = point(0);
        const double y = point(1);
        const double z = point(2);

        const double xComponent = std::pow(std::sin(10 * x), 2);
        const double yComponent = std::pow(std::sin(10 * y), 2);
        const double zComponent = std::pow(std::sin(10 * z), 2);

        return (xComponent + yComponent + zComponent) / 3;
    };

    grid.setDensitiesFunction(f);

    const auto file = std::filesystem::path(__FILE__).parent_path() / "output.msh";

    monad::saveGrid(grid, file.string(), true);

    std::cout << "Saved to " + file.string() << std::endl;

    return 0;
}
