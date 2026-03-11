#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/element/hex20.hpp"

namespace monad {

    namespace grid {

        /**
         * @brief Topology rules for 3D structured Hex20 grids.
         */
        struct Hex20Topology {
            using Element = fem::Hex20;
            using Resolution = std::array<std::size_t, Element::Dim>;
            using Size = std::array<double, Element::Dim>;
            using Point = typename Element::Point;
            using ElementList = std::array<std::size_t, Element::NumNodes>;

            /**
             * @brief Number of corner nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of corner nodes.
             */
            static std::size_t numCornerNodes(const Resolution& resolution) noexcept {
                return (resolution[0] + 1) * (resolution[1] + 1) * (resolution[2] + 1);
            }

            /**
             * @brief Number of x-edge midpoint nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of x-edge midpoint nodes.
             */
            static std::size_t numXMidNodes(const Resolution& resolution) noexcept {
                return resolution[0] * (resolution[1] + 1) * (resolution[2] + 1);
            }

            /**
             * @brief Number of y-edge midpoint nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of y-edge midpoint nodes.
             */
            static std::size_t numYMidNodes(const Resolution& resolution) noexcept {
                return (resolution[0] + 1) * resolution[1] * (resolution[2] + 1);
            }

            /**
             * @brief Number of z-edge midpoint nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of z-edge midpoint nodes.
             */
            static std::size_t numZMidNodes(const Resolution& resolution) noexcept {
                return (resolution[0] + 1) * (resolution[1] + 1) * resolution[2];
            }

            /**
             * @brief Number of nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of nodes.
             */
            static std::size_t numNodes(const Resolution& resolution) noexcept {
                return numCornerNodes(resolution) + numXMidNodes(resolution) + numYMidNodes(resolution) + numZMidNodes(resolution);
            }

            /**
             * @brief Number of periodic nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of periodic nodes.
             */
            static std::size_t numPeriodicNodes(const Resolution& resolution) noexcept {
                return 4 * resolution[0] * resolution[1] * resolution[2];
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
                const std::size_t nz = resolution[2];
                const std::size_t corners = numCornerNodes(resolution);
                const std::size_t xMids = numXMidNodes(resolution);
                const std::size_t yMids = numYMidNodes(resolution);

                const double dx = size[0] / static_cast<double>(nx);
                const double dy = size[1] / static_cast<double>(ny);
                const double dz = size[2] / static_cast<double>(nz);

                double x;
                double y;
                double z;

                if (index < corners) {
                    const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
                    const std::size_t indexInPlane = index % nodesPerPlane;
                    const std::size_t i = indexInPlane % (nx + 1);
                    const std::size_t j = indexInPlane / (nx + 1);
                    const std::size_t k = index / nodesPerPlane;

                    x = static_cast<double>(i) * dx;
                    y = static_cast<double>(j) * dy;
                    z = static_cast<double>(k) * dz;
                }
                else if (index < corners + xMids) {
                    index -= corners;

                    const std::size_t nodesPerPlane = nx * (ny + 1);
                    const std::size_t indexInPlane = index % nodesPerPlane;
                    const std::size_t i = indexInPlane % nx;
                    const std::size_t j = indexInPlane / nx;
                    const std::size_t k = index / nodesPerPlane;

                    x = (static_cast<double>(i) + 0.5) * dx;
                    y = static_cast<double>(j) * dy;
                    z = static_cast<double>(k) * dz;
                }
                else if (index < corners + xMids + yMids) {
                    index -= corners + xMids;

                    const std::size_t nodesPerPlane = (nx + 1) * ny;
                    const std::size_t indexInPlane = index % nodesPerPlane;
                    const std::size_t i = indexInPlane % (nx + 1);
                    const std::size_t j = indexInPlane / (nx + 1);
                    const std::size_t k = index / nodesPerPlane;

                    x = static_cast<double>(i) * dx;
                    y = (static_cast<double>(j) + 0.5) * dy;
                    z = static_cast<double>(k) * dz;
                }
                else {
                    index -= corners + xMids + yMids;

                    const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
                    const std::size_t indexInPlane = index % nodesPerPlane;
                    const std::size_t i = indexInPlane % (nx + 1);
                    const std::size_t j = indexInPlane / (nx + 1);
                    const std::size_t k = index / nodesPerPlane;

                    x = static_cast<double>(i) * dx;
                    y = static_cast<double>(j) * dy;
                    z = (static_cast<double>(k) + 0.5) * dz;
                }
                
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
                const std::size_t corners = numCornerNodes(resolution);
                const std::size_t xMids = numXMidNodes(resolution);
                const std::size_t yMids = numYMidNodes(resolution);

                const std::size_t i = index % nx;
                const std::size_t j = (index / nx) % ny;
                const std::size_t k = index / elementsPerPlane;

                const auto cornerNodeIndex = [nx, ny](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return kk * ((nx + 1) * (ny + 1)) + jj * (nx + 1) + ii;
                };

                const auto xMidNodeIndex = [nx, ny, corners](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return corners + kk * nx * (ny + 1) + jj * nx + ii;
                };

                const auto yMidNodeIndex = [nx, ny, corners, xMids](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return corners + xMids + kk * (nx + 1) * ny + jj * (nx + 1) + ii;
                };

                const auto zMidNodeIndex = [nx, ny, corners, xMids, yMids](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return corners + xMids + yMids + kk * (nx + 1) * (ny + 1) + jj * (nx + 1) + ii;
                };

                return {
                    cornerNodeIndex(i, j, k),
                    cornerNodeIndex(i + 1, j, k),
                    cornerNodeIndex(i + 1, j + 1, k),
                    cornerNodeIndex(i, j + 1, k),
                    cornerNodeIndex(i, j, k + 1),
                    cornerNodeIndex(i + 1, j, k + 1),
                    cornerNodeIndex(i + 1, j + 1, k + 1),
                    cornerNodeIndex(i, j + 1, k + 1),
                    xMidNodeIndex(i, j, k),
                    yMidNodeIndex(i + 1, j, k),
                    xMidNodeIndex(i, j + 1, k),
                    yMidNodeIndex(i, j, k),
                    xMidNodeIndex(i, j, k + 1),
                    yMidNodeIndex(i + 1, j, k + 1),
                    xMidNodeIndex(i, j + 1, k + 1),
                    yMidNodeIndex(i, j, k + 1),
                    zMidNodeIndex(i, j, k),
                    zMidNodeIndex(i + 1, j, k),
                    zMidNodeIndex(i + 1, j + 1, k),
                    zMidNodeIndex(i, j + 1, k)
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
                const std::size_t numElements = nx * ny * nz;
                const std::size_t elementsPerPlane = nx * ny;

                const std::size_t i = index % nx;
                const std::size_t j = (index / nx) % ny;
                const std::size_t k = index / elementsPerPlane;

                const auto cornerNodeIndex = [nx, ny, nz](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return (kk % nz) * (nx * ny) + (jj % ny) * nx + (ii % nx);
                };

                const auto xMidNodeIndex = [nx, ny, nz, numElements](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return numElements + (kk % nz) * (nx * ny) + (jj % ny) * nx + (ii % nx);
                };

                const auto yMidNodeIndex = [nx, ny, nz, numElements](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return 2 * numElements + (kk % nz) * (nx * ny) + (jj % ny) * nx + (ii % nx);
                };

                const auto zMidNodeIndex = [nx, ny, nz, numElements](std::size_t ii, std::size_t jj, std::size_t kk) {
                    return 3 * numElements + (kk % nz) * (nx * ny) + (jj % ny) * nx + (ii % nx);
                };

                return {
                    cornerNodeIndex(i, j, k),
                    cornerNodeIndex(i + 1, j, k),
                    cornerNodeIndex(i + 1, j + 1, k),
                    cornerNodeIndex(i, j + 1, k),
                    cornerNodeIndex(i, j, k + 1),
                    cornerNodeIndex(i + 1, j, k + 1),
                    cornerNodeIndex(i + 1, j + 1, k + 1),
                    cornerNodeIndex(i, j + 1, k + 1),
                    xMidNodeIndex(i, j, k),
                    yMidNodeIndex(i + 1, j, k),
                    xMidNodeIndex(i, j + 1, k),
                    yMidNodeIndex(i, j, k),
                    xMidNodeIndex(i, j, k + 1),
                    yMidNodeIndex(i + 1, j, k + 1),
                    xMidNodeIndex(i, j + 1, k + 1),
                    yMidNodeIndex(i, j, k + 1),
                    zMidNodeIndex(i, j, k),
                    zMidNodeIndex(i + 1, j, k),
                    zMidNodeIndex(i + 1, j + 1, k),
                    zMidNodeIndex(i, j + 1, k)
                };
            }
        };

    } // namespace grid

} // namespace monad
