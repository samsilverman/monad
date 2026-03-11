#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/element/hex8.hpp"

namespace monad {

    namespace grid {

        /**
         * @brief Topology rules for 3D structured Hex8 grids.
         */
        struct Hex8Topology {
            using Element = fem::Hex8;
            using Resolution = std::array<std::size_t, Element::Dim>;
            using Size = std::array<double, Element::Dim>;
            using Point = typename Element::Point;
            using ElementList = std::array<std::size_t, Element::NumNodes>;

            /// @brief Number of nodes.
            static std::size_t numNodes(const Resolution& resolution) noexcept {
                return (resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1);
            }

            /// @brief Number of periodic nodes.
            static std::size_t numPeriodicNodes(const Resolution& resolution) noexcept {
                return resolution[0] * resolution[1] * resolution[2];
            }

            /**
             * @brief Coordinates for a specific node.
             *
             * @param[in] index Node index.
             * @param[in] resolution How many cells in each dimension.
             * @param[in] size Physical lengths in each dimension.
             *
             * @returns Coordinates for a specific node.
             */
            static Point node(std::size_t index, const Resolution& resolution, const Size& size) noexcept {
                const std::size_t nx = resolution[0];
                const std::size_t ny = resolution[1];
                const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
                const std::size_t indexInPlane = index % nodesPerPlane;

                const std::size_t i = indexInPlane % (nx + 1);
                const std::size_t j = indexInPlane / (nx + 1);
                const std::size_t k = index / nodesPerPlane;

                const double dx = size[0] / static_cast<double>(nx);
                const double dy = size[1] / static_cast<double>(ny);
                const double dz = size[2] / static_cast<double>(resolution[2]);

                const double x = static_cast<double>(i) * dx;
                const double y = static_cast<double>(j) * dy;
                const double z = static_cast<double>(k) * dz;

                return Point(x, y, z);
            }

            /**
             * @brief Node indices for a specific element.
             *
             * @param[in] index Element index.
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Node indices for a specific element.
             */
            static ElementList element(std::size_t index, const Resolution& resolution) noexcept {
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

            /**
             * @brief Periodic node indices for a specific element.
             *
             * @param[in] index Element index.
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Periodic node indices for a specific element.
             */
            static ElementList periodicElement(std::size_t index, const Resolution& resolution) noexcept {
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
        };

    } // namespace grid

} // namespace monad
