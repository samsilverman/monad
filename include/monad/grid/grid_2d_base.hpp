#pragma once

#include <fstream>
#include <stdexcept>
#include <cstddef>
#include <vector>
#include <string>
#include <sstream>
#include "monad/grid/grid_base.hpp"

namespace monad {

    /**
     * @brief Curiously recurring template pattern (CRTP) base class for 2D grid meshes.
     *
     * @tparam Derived Concrete grid mesh class (CRTP).
     * @tparam Element Concrete element class (e.g. Quad4).
     *
     * @note `Derived` classes must provide:
     *
     * - `std::size_t numNodes() const noexcept`
     *
     * - `std::size_t numPeriodicNodes() const noexcept`
     *
     * - `Point node(std::size_t index) const`
     *
     * - `ElementList element(std::size_t index) const`
     *
     * - `ElementList periodicElement(std::size_t index) const`
     */
    template <class Derived, class Element>
    class Grid2dBase : public GridBase<Derived, Element> {
    public:
        static_assert(Element::Dim == 2, "Element spatial dimension must be 2.");

        using GridBase<Derived, Element>::GridBase;
        using Resolution = typename GridBase<Derived, Element>::Resolution;
        using DensityList = typename GridBase<Derived, Element>::DensityList;

        /**
         * @brief Translates the element densities periodically.
         *
         * @param[in] shift Shift in each directions.
         */
        void translate(const Resolution &shift) noexcept {
            DensityList shiftedDensities(this->numElements());

            const std::size_t nx = this->resolution_[0];
            const std::size_t ny = this->resolution_[1];

            for (std::size_t i = 0; i < nx; ++i) {
                for (std::size_t j = 0; j < ny; ++j) {
                    const std::size_t oldIndex = j * nx + i;

                    const std::size_t iNew = (i + shift[0]) % nx;
                    const std::size_t jNew = (j + shift[1]) % ny;

                    std::size_t newIndex = jNew * nx + iNew;

                    shiftedDensities[newIndex] = this->densities_[oldIndex];
                }
            }

            this->setDensities(shiftedDensities);
        }

        /// @brief Grid area.
        double area() const noexcept {
            return this->measure();
        }

        /**
         * @brief Set the material densities using a CSV file.
         *
         * @param[in] file CSV file.
         *
         * @throws std::runtime_error if `file` cannot be opened.
         * @throws std::runtime_error if `file` contains no data.
         * @throws std::runtime_error if `file` contains non-numeric data.
         * @throws std::runtime_error if `file` contains numbers outside the range [0, 1].
         * @throws std::runtime_error if `file` contains rows/columns not matching the grid resolution.
         *
         * @note The CSV data is interpreted as a 2D array with the origin
         * located at the bottom-left corner.
         */
        void setDensitiesFile(const std::string &file) {
            std::ifstream ifs(file, std::ios::in | std::ios::binary);
            if (!ifs.is_open()) {
                throw std::runtime_error("Could not open " + file + " for reading.");
            }

            const std::size_t nx = this->resolution_[0];
            const std::size_t ny = this->resolution_[1];

            std::vector<std::vector<double>> rows;

            std::string line;

            // Load csv data
            while (std::getline(ifs, line)) {
                if (line.empty()) {
                    continue;
                }

                std::istringstream iss(line);
                std::string cell;

                std::vector<double> row;
                row.reserve(nx);

                while (std::getline(iss, cell, ',')) {
                    if (cell.empty()) {
                        continue;
                    }

                    double value;
                    try {
                        value = std::stod(cell);
                    }
                    catch (const std::invalid_argument &) {
                        throw std::runtime_error("File " + file + " contains non-numeric data.");
                    }
                    if (value < 0.0 || value > 1.0) {
                        throw std::runtime_error("File " + file + " contains data outside the range [0,1].");
                    }
                    row.push_back(value);
                }

                if (row.empty()) {
                    continue;
                }

                if (row.size() != nx) {
                    throw std::runtime_error("File " + file + " number of columns (" + std::to_string(row.size()) + ") does not equal grid x-resolution (" + std::to_string(nx) + ").");
                }

                rows.push_back(row);
            }

            if (rows.size() != ny) {
                throw std::runtime_error("File " + file + " number of rows (" + std::to_string(rows.size()) + ") does not equal grid y-resolution (" + std::to_string(ny) + ").");
            }

            // Flatten csv data
            for (std::size_t i = 0; i < ny; ++i) {
                for (std::size_t j = 0; j < nx; ++j) {
                    std::size_t rowStart = nx * (ny - 1 - i);
                    this->setDensity(rowStart + j, rows[i][j]);
                }
            }
        }
    };

} // namespace monad
