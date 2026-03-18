#include "monad/grid/quad8_topology.hpp"

namespace monad {

    namespace grid {

        std::size_t Quad8Topology::numCornerNodes(const Resolution& resolution) noexcept {
            return (resolution[0] + 1) * (resolution[1] + 1);
        }

        std::size_t Quad8Topology::numXMidNodes(const Resolution& resolution) noexcept {
            return resolution[0] * (resolution[1] + 1);
        }

        std::size_t Quad8Topology::numYMidNodes(const Resolution& resolution) noexcept {
            return (resolution[0] + 1) * resolution[1];
        }

        std::size_t Quad8Topology::numNodes(const Resolution& resolution) noexcept {
            return numCornerNodes(resolution) + numXMidNodes(resolution) + numYMidNodes(resolution);
        }

        std::size_t Quad8Topology::numPeriodicNodes(const Resolution& resolution) noexcept {
            return 3 * resolution[0] * resolution[1];
        }

        Quad8Topology::Point Quad8Topology::node(std::size_t index, const Resolution& resolution, const Size& size) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t corners = numCornerNodes(resolution);
            const std::size_t xMids = numXMidNodes(resolution);

            const double dx = size[0] / static_cast<double>(nx);
            const double dy = size[1] / static_cast<double>(ny);

            double x;
            double y;

            if (index < corners) {
                const std::size_t i = index % (nx + 1);
                const std::size_t j = index / (nx + 1);

                x = static_cast<double>(i) * dx;
                y = static_cast<double>(j) * dy;
            }
            else if (index < corners + xMids) {
                index -= corners;

                const std::size_t i = index % nx;
                const std::size_t j = index / nx;

                x = (static_cast<double>(i) + 0.5) * dx;
                y = static_cast<double>(j) * dy;
            }
            else {
                index -= corners + xMids;

                const std::size_t i = index % (nx + 1);
                const std::size_t j = index / (nx + 1);

                x = static_cast<double>(i) * dx;
                y = (static_cast<double>(j) + 0.5) * dy;
            }

            return Point(x, y);
        }

        Quad8Topology::ElementList Quad8Topology::element(std::size_t index, const Resolution& resolution) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t corners = numCornerNodes(resolution);
            const std::size_t xMids = numXMidNodes(resolution);

            const std::size_t i = index % nx;
            const std::size_t j = index / nx;

            const auto cornerNodeIndex = [nx](std::size_t ii, std::size_t jj) {
                return jj * (nx + 1) + ii;
            };

            const auto xMidNodeIndex = [nx, corners](std::size_t ii, std::size_t jj) {
                return corners + jj * nx + ii;
            };

            const auto yMidNodeIndex = [nx, corners, xMids](std::size_t ii, std::size_t jj) {
                return corners + xMids + jj * (nx + 1) + ii;
            };

            return {
                cornerNodeIndex(i, j),
                cornerNodeIndex(i + 1, j),
                cornerNodeIndex(i + 1, j + 1),
                cornerNodeIndex(i, j + 1),
                xMidNodeIndex(i, j),
                yMidNodeIndex(i + 1, j),
                xMidNodeIndex(i, j + 1),
                yMidNodeIndex(i, j)
            };
        }

        Quad8Topology::ElementList Quad8Topology::periodicElement(std::size_t index, const Resolution& resolution) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t numElements = nx * ny;

            const std::size_t i = index % nx;
            const std::size_t j = index / nx;

            const auto cornerNodeIndex = [nx, ny](std::size_t ii, std::size_t jj) {
                return (jj % ny) * nx + (ii % nx);
            };

            const auto xMidNodeIndex = [nx, ny, numElements](std::size_t ii, std::size_t jj) {
                return numElements + (jj % ny) * nx + (ii % nx);
            };

            const auto yMidNodeIndex = [nx, ny, numElements](std::size_t ii, std::size_t jj) {
                return 2 * numElements + (jj % ny) * nx + (ii % nx);
            };

            return {
                cornerNodeIndex(i, j),
                cornerNodeIndex(i + 1, j),
                cornerNodeIndex(i + 1, j + 1),
                cornerNodeIndex(i, j + 1),
                xMidNodeIndex(i, j),
                yMidNodeIndex(i + 1, j),
                xMidNodeIndex(i, j + 1),
                yMidNodeIndex(i, j)
            };
        }

    } // namespace grid

} // namespace monad
