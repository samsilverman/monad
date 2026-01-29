/**
 * @brief Command-line tool for generating a 2D mesh.
 *
 * This program creates a 2D Quad8 grid mesh and
 * saves it to a Gmsh-compatible file.
 *
 * The Gmsh file is written to:
 *
 *      `/path/to/apps/1_2DGrid/output.gmsh`
 *
 * @param nx Number of grid cells in the x-direction.
 * @param ny Number of grid cells in the y-direction.
 * @param lx Length of grid in the x-direction (default=1).
 * @param ly Length of grid in the y-direction (default=1).
 *
 * Usage:
 *      $ 1_2DGrid <nx> <ny> [lx] [ly]
 *
 * Example:
 *      $ 1_2DGrid 10 5 1.0 0.5
 *      Saved to /path/to/apps/1_2DGrid/output.gmsh
 */
#include <iostream>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include "monad/monad.hpp"

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <nx> <ny> [lx] [ly]\n";
        return 1;
    }

    std::size_t nx;
    std::size_t ny;
    double lx = 1.0;
    double ly = 1.0;

    try {
        nx = std::stoul(argv[1]);
        ny = std::stoul(argv[2]);

        if (argc > 3) {
            lx = std::stod(argv[3]);
        }
        if (argc > 4) {
            ly = std::stod(argv[4]);
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "1_2DGrid: Invalid argument: " << e.what() << std::endl;
        return 1;
    }

    if (nx == 0) {
        std::cerr << "1_2DGrid: nx must be positive.\n";
        return 1;
    }
    if (ny == 0) {
        std::cerr << "1_2DGrid: ny must be positive.\n";
        return 1;
    }
    if (lx <= 0) {
        std::cerr << "1_2DGrid: lx (" << lx <<  ") must be positive.\n";
        return 1;
    }
    if (ly <= 0) {
        std::cerr << "1_2DGrid: ly (" << ly <<  ") must be positive.\n";
        return 1;
    }

    const monad::Quad8Grid grid({nx, ny}, {lx, ly});

    const auto file = std::filesystem::path(__FILE__).parent_path() / "output.msh";

    monad::saveGrid(grid, file.string());

    std::cout << "Saved to " + file.string() << std::endl;

    return 0;
}
