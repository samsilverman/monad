#pragma once

#include <cstddef>
#include <ostream>

namespace monad {

    namespace io {

        namespace gmsh {

            /**
             * @brief Writes the Gmsh file `$Elements` section.
             *
             * @tparam Grid Grid class (e.g. Quad4Grid).
             *
             * @param[in,out] os Output stream.
             * @param[in] grid Periodic unit cell grid.
             *
             * @note Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
             */
            template <typename Grid>
            void writeGmshElements(std::ostream &os, const Grid &grid) noexcept {
                static_assert(Grid::Dim == 2 || Grid::Dim == 3, "Grid spatial dimension must be 2 or 3.");

                const std::size_t numEntityBlocks = 1;
                const std::size_t numElements = grid.numElements();
                const std::size_t minElementTag = 1;
                const std::size_t maxElementTag = numElements;
                constexpr int entityDim = Grid::Dim;
                const int entityTag = 1;
                const int elementType = Grid::Element::gmshElementType();
                const std::size_t numElementsInBlock = numElements;

                // Section header
                os << "$Elements\n"
                   << numEntityBlocks << " " << numElements << " " << minElementTag << " " << maxElementTag << "\n"
                   << entityDim << " " << entityTag << " " << elementType << " " << numElementsInBlock << "\n";

                // Section body
                const auto ordering = Grid::Element::gmshNodeOrdering();

                for (std::size_t i = 0; i < numElements; ++i) {
                    const std::size_t elementTag = i + 1;

                    os << elementTag;

                    const auto element = grid.element(i);

                    for (std::size_t j : ordering) {
                        const std::size_t nodeTag = element[j] + 1;

                        os << " " << nodeTag;
                    }

                    os << "\n";
                }

                // Section footer
                os << "$EndElements";
            }

        } // namespace gmsh

    } // namespace io

} // namespace monad
