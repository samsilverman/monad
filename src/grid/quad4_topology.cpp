#include "monad/grid/quad4_topology.hpp"

namespace monad {

    namespace grid {

        std::size_t Quad4Topology::numNodes(const Resolution& resolution) noexcept {
            return (resolution[0] + 1) * (resolution[1] + 1);
        }

        std::size_t Quad4Topology::numPeriodicNodes(const Resolution& resolution) noexcept {
            return resolution[0] * resolution[1];
        }

        Quad4Topology::Point Quad4Topology::node(std::size_t index, const Resolution& resolution, const Size& size) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t i = index % (nx + 1);
            const std::size_t j = index / (nx + 1);

            const double dx = size[0] / static_cast<double>(nx);
            const double dy = size[1] / static_cast<double>(ny);

            const double x = static_cast<double>(i) * dx;
            const double y = static_cast<double>(j) * dy;

            return Point(x, y);
        }

        Quad4Topology::ElementList Quad4Topology::element(std::size_t index, const Resolution& resolution) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t i = index % nx;
            const std::size_t j = index / nx;

            const auto nodeIndex = [nx](std::size_t ii, std::size_t jj) {
                return jj * (nx + 1) + ii;
            };

            return {
                nodeIndex(i, j),
                nodeIndex(i + 1, j),
                nodeIndex(i + 1, j + 1),
                nodeIndex(i, j + 1)
            };
        }

        Quad4Topology::ElementList Quad4Topology::periodicElement(std::size_t index, const Resolution& resolution) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t i = index % nx;
            const std::size_t j = index / nx;

            const auto nodeIndex = [nx, ny](std::size_t ii, std::size_t jj) {
                return (jj % ny) * nx + (ii % nx);
            };

            return {
                nodeIndex(i, j),
                nodeIndex(i + 1, j),
                nodeIndex(i + 1, j + 1),
                nodeIndex(i, j + 1)
            };
        }

    } // namespace grid

} // namespace monad
