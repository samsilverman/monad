#pragma once

#include <string>
#include <filesystem>
#include <stdexcept>
#include <fstream>
#include <Eigen/Core>
#include "monad/io/gmsh/write_gmsh_header.hpp"
#include "monad/io/gmsh/write_gmsh_nodes.hpp"
#include "monad/io/gmsh/write_gmsh_elements.hpp"
#include "monad/io/gmsh/write_gmsh_nodal_field.hpp"

namespace monad {
    /**
     * @brief Writes a grid and nodal field to a Gmsh file.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Derived Eigen matrix type.
     *
     * @param[in] grid Periodic unit cell grid.
     * @param[in] file File.
     * @param[in] field An Nx1 (scalar), Nx2 or Nx3 (vector) nodal field.
     * @param[in] name Optional name tag for the nodal field.
     *
     * @throws std::invalid_argument if the `file` extension is not `.msh`.
     * @throws std::invalid_argument if the `field` size (N) does not equal the number of `grid` nodes.
     * @throws std::runtime_error if the `file` cannot be opened for writing.
     */
    template <typename Grid, typename Derived>
    void saveGridAndField(const Grid &grid, const Eigen::MatrixBase<Derived>& field, const std::string &file, const std::string &name = "") {
        static_assert(Grid::Dim == 2 || Grid::Dim == 3, "Grid spatial dimension must be 2 or 3.");

        if (std::filesystem::path(file).extension() != ".msh") {
            throw std::invalid_argument("File extension must be \".msh\".");
        }

        if (field.rows() != static_cast<int>(grid.numNodes())) {
            throw std::invalid_argument("Field size (" + std::to_string(field.rows()) + ") must equal number of grid nodes (" + std::to_string(grid.numNodes()) + ").");
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
        ofs << "\n\n";
        io::gmsh::writeGmshNodalField(ofs, field, name);
        ofs << "\n";
    }

} // namespace monad
