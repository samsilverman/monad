#pragma once

#include <string>
#include "monad/field/field_aliases.hpp" 
    
namespace monad {

    /**
     * @brief Constructs a 2D density field from a CSV file.
     *
     * @param[in] file Path to the CSV file.
     *
     * @returns 2D density field defined by the `file` entries.
     *
     * @throws std::runtime_error if `file` cannot be opened.
     * @throws std::runtime_error if `file` is empty.
     * @throws std::runtime_error if `file` contains non-numeric values.
     * @throws std::runtime_error if `file` contains any value outside the range [0,1].
     * @throws std::runtime_error if `file` contains rows with inconsistent lengths.
     *
     * @note CSV entries are interpreted as a 2D grid with the origin at the
     * bottom-left corner.
     */
    DensityField2d makeDensityFieldFromCsv(const std::string& file);

} // namespace monad
