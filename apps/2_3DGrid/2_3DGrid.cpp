/**
 * @brief Command-line tool for generating a 3D mesh.
 *
 * This program creates a 3D Hex8 grid mesh with random
 * densities and saves it to a Gmsh-compatible file.
 *
 * The Gmsh file is written to:
 *
 *      `/path/to/apps/2_3DGrid/output.gmsh`
 *
 * @param nx Number of grid cells in the x-direction.
 * @param ny Number of grid cells in the y-direction.
 * @param nz Number of grid cells in the z-direction.
 * @param lx Length of grid in the x-direction (default=1).
 * @param ly Length of grid in the y-direction (default=1).
 * @param lz Length of grid in the z-direction (default=1).
 * @param seed Seed for random number generation (default=1234).
 * Usage:
 *      $ 2_3DGrid <nx> <ny> <nz> [lx] [ly] [lz] [seed]
 *
 * Example:
 *      $ 2_3DGrid 15 10 5 0.3 0.2 0.1 1234
 *      Saved to /path/to/apps/2_3DGrid/output.gmsh
 */
#include <iostream>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include "monad/monad.hpp"

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 8) {
        std::cerr << "Usage: " << argv[0] << " <nx> <ny> <nz> [lx] [ly] [lz] [seed]\n";
        return 1;
    }

    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
    double lx = 1.0;
    double ly = 1.0;
    double lz = 1.0;
    int seed = 1234;

    try {
        nx = std::stoul(argv[1]);
        ny = std::stoul(argv[2]);
        nz = std::stoul(argv[3]);

        if (argc > 4) {
            lx = std::stod(argv[4]);
        }
        if (argc > 5) {
            ly = std::stod(argv[5]);
        }
        if (argc > 6) {
            lz = std::stod(argv[6]);
        }
        if (argc > 7) {
            seed = std::stoi(argv[7]);
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "2_3DGrid: Invalid argument: " << e.what() << std::endl;
        return 1;
    }

    if (nx <= 0) {
        std::cerr << "2_3DGrid: nx must be positive.\n";
        return 1;
    }
    if (ny <= 0) {
        std::cerr << "2_3DGrid: ny must be positive.\n";
        return 1;
    }
    if (nz <= 0) {
        std::cerr << "2_3DGrid: nz must be positive.\n";
        return 1;
    }
    if (lx <= 0) {
        std::cerr << "2_3DGrid: lx (" << lx <<  ") must be positive.\n";
        return 1;
    }
    if (ly <= 0) {
        std::cerr << "2_3DGrid: ly (" << ly <<  ") must be positive.\n";
        return 1;
    }
    if (lz <= 0) {
        std::cerr << "2_3DGrid: ly (" << lz <<  ") must be positive.\n";
        return 1;
    }
    if (seed <= 0) {
        std::cerr << "2_3DGrid: seed (" << seed <<  ") must be positive.\n";
        return 1;
    }

    monad::Hex8Grid grid({nx, ny, nz}, {lx, ly, lz});
    grid.setDensitiesRandom(seed);

    const auto file = std::filesystem::path(__FILE__).parent_path() / "output.msh";

    monad::saveGrid(grid, file.string(), true);

    std::cout << "Saved to " + file.string() << std::endl;

    return 0;
}
