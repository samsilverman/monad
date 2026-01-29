#pragma once

#include <cstddef>
#include <string>
#include <ostream>
#include "monad/detail/constants.hpp"

namespace monad {

    namespace io {

        namespace gmsh {

            /**
             * @brief Writes the Gmsh file `$ElementData` section containing per-element material densities.
             *
             * @tparam Grid Grid class (e.g. Quad4Grid).
             *
             * @param[in,out] os Output stream.
             * @param[in] grid Periodic unit cell grid.
             *
             * @note Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
             */
            template <typename Grid>
            void writeGmshDensities(std::ostream &os, const Grid &grid) noexcept {
                const int numStringTags = 1;
                const std::string stringTag = "Density";
                const int numRealTags = 0;
                // 3 tags: timestep, data dimension, number of entries.
                const int numIntegerTags = 3;

                // Section header
                os << "$ElementData\n"
                   << numStringTags << "\n"
                   << "\"" << stringTag << "\"\n"
                   << numRealTags << "\n"
                   << numIntegerTags << "\n"
                   << "0\n"
                   << "1\n"
                   << static_cast<int>(grid.numElements()) << "\n";

                // Section body
                for (std::size_t i = 0; i < grid.numElements(); ++i) {
                    const std::size_t elementTag = i + 1;

                    os << elementTag;

                    double value = grid.getDensity(i);

                    // Report 0 rather than NUMERICAL_ZERO for cleaner output
                    if (value <= NUMERICAL_ZERO) {
                        value = 0.0;
                    }

                    os << " " << value << "\n";
                }

                // Section footer
                os << "$EndElementData";
            }

        } // namespace gmsh

    } // namespace io

} // namespace monad
