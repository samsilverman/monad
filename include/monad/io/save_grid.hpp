#pragma once

#include <string>
#include <filesystem>
#include <stdexcept>
#include <fstream>
#include "monad/io/gmsh/write_gmsh_header.hpp"
#include "monad/io/gmsh/write_gmsh_nodes.hpp"
#include "monad/io/gmsh/write_gmsh_elements.hpp"
#include "monad/io/gmsh/write_gmsh_densities.hpp"

namespace monad {
    /**
     * @brief Writes a grid to a Gmsh file.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     *
     * @param[in] grid Grid.
     * @param[in] file Path to the Gmsh file.
     *
     * @throws std::invalid_argument if `file` does not have the `.msh` extension.
     * @throws std::runtime_error if the `file` cannot be opened for writing.
     */
    template <typename Grid>
    void saveGrid(const Grid &grid, const std::string &file) {
        static_assert(Grid::Dim == 2 || Grid::Dim == 3, "Grid spatial dimension must be 2 or 3.");

        if (std::filesystem::path(file).extension() != ".msh") {
            throw std::invalid_argument("File extension must be \".msh\".");
        }

        std::ofstream ofs(file, std::ios::trunc);
        if (!ofs.is_open()) {
            throw std::runtime_error("Could not open " + file + " for writing.");
        }

        io::gmsh::writeGmshHeader(ofs);
        ofs << "\n\n";
        io::gmsh::writeGmshNodes(ofs, grid);
        ofs << "\n\n";
        io::gmsh::writeGmshElements(ofs, grid);
        ofs << "\n";
    }

} // namespace monad
