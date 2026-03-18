#pragma once

#include <cstddef>
#include <string>
#include <ostream>
#include "monad/field/density_field.hpp"
#include "monad/detail/constants.hpp"

namespace monad {

    namespace io {

        namespace gmsh {

            /**
             * @brief Writes the Gmsh `$ElementData` section for per-element material densities.
             *
             * @tparam D Spatial dimension (2 or 3).
             *
             * @param[in,out] os Output stream.
             * @param[in] densityField Per-element density field.
             *
             * @note Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
             */
            template <int D>
            void writeGmshDensities(std::ostream &os, const field::DensityField<D> &densityField) noexcept {
                static_assert(D == 2 || D == 3, "Density field spatial dimension must be 2 or 3.");

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
                   << static_cast<int>(densityField.numElements()) << "\n";

                // Section body
                for (std::size_t i = 0; i < densityField.numElements(); ++i) {
                    const std::size_t elementTag = i + 1;

                    os << elementTag;

                    double value = densityField.getDensity(i);

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
