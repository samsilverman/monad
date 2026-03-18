#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "monad/field/make_density_field_from_csv.hpp"

namespace monad {

    DensityField2d makeDensityFieldFromCsv(const std::string& file) {
        std::ifstream ifs(file, std::ios::in | std::ios::binary);
        if (!ifs.is_open()) {
            throw std::runtime_error("Could not open " + file + " for reading.");
        }

        std::size_t nx;
        std::vector<std::vector<double>> rows;
        std::string line;

        while (std::getline(ifs, line)) {
            if (line.empty()) {
                continue;
            }

            std::istringstream iss(line);
            std::string cell;
            std::vector<double> row;

            while (std::getline(iss, cell, ',')) {
                if (cell.empty()) {
                    continue;
                }

                double value = 0.0;
                try {
                    value = std::stod(cell);
                }
                catch (const std::invalid_argument&) {
                    throw std::runtime_error("File " + file + " contains non-numeric values.");
                }

                if (value < 0.0 || value > 1.0) {
                    throw std::runtime_error("File " + file + " contains values outside the range [0,1].");
                }

                row.push_back(value);
            }

            if (row.empty()) {
                continue;
            }

            if (rows.size() == 0) {
                nx = row.size();
            }

            if (row.size() != nx) {
                throw std::runtime_error("File " + file + " number of columns (" + std::to_string(row.size()) + ") does not equal density field x-resolution (" + std::to_string(nx) + ").");
            }

            rows.push_back(row);
        }

        const std::size_t ny = rows.size();

        if (ny == 0) {
            throw std::runtime_error("File " + file + " is empty.");
        }

        DensityField2d densityField({nx, ny});

        for (std::size_t i = 0; i < ny; ++i) {
            for (std::size_t j = 0; j < nx; ++j) {
                const std::size_t rowStart = nx * (ny - 1 - i);
                densityField.setDensity(rowStart + j, rows[i][j]);
            }
        }

        return densityField;
    }

} // namespace monad
