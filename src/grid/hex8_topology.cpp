#include "monad/grid/hex8_topology.hpp"

namespace monad {

    namespace grid {

        std::size_t Hex8Topology::numNodes(const Resolution& resolution) noexcept {
            return (resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1);
        }

        std::size_t Hex8Topology::numPeriodicNodes(const Resolution& resolution) noexcept {
            return resolution[0] * resolution[1] * resolution[2];
        }

        Hex8Topology::Point Hex8Topology::node(std::size_t index, const Resolution& resolution, const Size& size) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t nz = resolution[2];
            const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
            const std::size_t indexInPlane = index % nodesPerPlane;

            const std::size_t i = indexInPlane % (nx + 1);
            const std::size_t j = indexInPlane / (nx + 1);
            const std::size_t k = index / nodesPerPlane;

            const double dx = size[0] / static_cast<double>(nx);
            const double dy = size[1] / static_cast<double>(ny);
            const double dz = size[2] / static_cast<double>(nz);

            const double x = static_cast<double>(i) * dx;
            const double y = static_cast<double>(j) * dy;
            const double z = static_cast<double>(k) * dz;

            return Point(x, y, z);
        }

        Hex8Topology::ElementList Hex8Topology::element(std::size_t index, const Resolution& resolution) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t elementsPerPlane = nx * ny;

            const std::size_t i = index % nx;
            const std::size_t j = (index / nx) % ny;
            const std::size_t k = index / elementsPerPlane;

            const auto nodeIndex = [nx, ny](std::size_t ii, std::size_t jj, std::size_t kk) {
                return kk * ((nx + 1) * (ny + 1)) + jj * (nx + 1) + ii;
            };

            return {
                nodeIndex(i, j, k),
                nodeIndex(i + 1, j, k),
                nodeIndex(i + 1, j + 1, k),
                nodeIndex(i, j + 1, k),
                nodeIndex(i, j, k + 1),
                nodeIndex(i + 1, j, k + 1),
                nodeIndex(i + 1, j + 1, k + 1),
                nodeIndex(i, j + 1, k + 1)
            };
        }

        Hex8Topology::ElementList Hex8Topology::periodicElement(std::size_t index, const Resolution& resolution) noexcept {
            const std::size_t nx = resolution[0];
            const std::size_t ny = resolution[1];
            const std::size_t nz = resolution[2];
            const std::size_t elementsPerPlane = nx * ny;

            const std::size_t i = index % nx;
            const std::size_t j = (index / nx) % ny;
            const std::size_t k = index / elementsPerPlane;

            const auto nodeIndex = [nx, ny, nz](std::size_t ii, std::size_t jj, std::size_t kk) {
                return (kk % nz) * (nx * ny) + (jj % ny) * nx + (ii % nx);
            };

            return {
                nodeIndex(i, j, k),
                nodeIndex(i + 1, j, k),
                nodeIndex(i + 1, j + 1, k),
                nodeIndex(i, j + 1, k),
                nodeIndex(i, j, k + 1),
                nodeIndex(i + 1, j, k + 1),
                nodeIndex(i + 1, j + 1, k + 1),
                nodeIndex(i, j + 1, k + 1)
            };
        }

    } // namespace grid

} // namespace monad
