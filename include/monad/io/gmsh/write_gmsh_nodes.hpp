#pragma once

#include <cstddef>
#include <ostream>

namespace monad {

    namespace io {

        namespace gmsh {

            /**
             * @brief Writes the Gmsh file `$Nodes` section.
             *
             * @tparam Grid Grid class (e.g. Quad4Grid).
             *
             * @param[in,out] os Output stream.
             * @param[in] grid Periodic unit cell grid.
             *
             * @note Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
             */
            template <typename Grid>
            void writeGmshNodes(std::ostream &os, const Grid &grid) noexcept {
                static_assert(Grid::Dim == 2 || Grid::Dim == 3, "Grid spatial dimension must be 2 or 3.");

                const std::size_t numEntityBlocks = 1;
                const std::size_t numNodes = grid.numNodes();
                const std::size_t minNodeTag = 1;
                const std::size_t maxNodeTag = numNodes;
                constexpr int entityDim = Grid::Dim;
                const int entityTag = 1;
                const int parametric = 0;
                const std::size_t numNodesInBlock = numNodes;

                // Section header
                os << "$Nodes\n"
                   << numEntityBlocks << " " << numNodes << " " << minNodeTag << " " << maxNodeTag << "\n"
                   << entityDim << " " << entityTag << " " << parametric << " " << numNodesInBlock << "\n";

                for (std::size_t nodeTag = 1; nodeTag <= numNodes; ++nodeTag) {
                    os << nodeTag << "\n";
                }

                // Section body
                for (std::size_t i = 0; i < numNodes; ++i) {
                    const auto node = grid.node(i);

                    const double x = node(0);
                    const double y = node(1);
                    double z = 0.0;

                    if constexpr (entityDim == 3) {
                        z = node(2);
                    }

                    os << x << " " << y << " " << z << "\n";
                }

                // Section footer
                os << "$EndNodes";
            }

        } // namespace gmsh

    } // namespace io

} // namespace monad
