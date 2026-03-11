#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/element/quad8.hpp"

namespace monad {

    namespace grid {

        /**
         * @brief Topology rules for 2D structured Quad8 grids.
         */
        struct Quad8Topology {
            using Element = fem::Quad8;
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
                return (resolution[0] + 1) * (resolution[1] + 1);
            }

            /**
             * @brief Number of horizontal edge midpoint nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of horizontal edge midpoint nodes.
             */
            static std::size_t numXMidNodes(const Resolution& resolution) noexcept {
                return resolution[0] * (resolution[1] + 1);
            }

            /**
             * @brief Number of vertical edge midpoint nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of vertical edge midpoint nodes.
             */
            static std::size_t numYMidNodes(const Resolution& resolution) noexcept {
                return (resolution[0] + 1) * resolution[1];
            }

            /**
             * @brief Number of nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of nodes.
             */
            static std::size_t numNodes(const Resolution& resolution) noexcept {
                return numCornerNodes(resolution) + numXMidNodes(resolution) + numYMidNodes(resolution);
            }

            /**
             * @brief Number of periodic nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of periodic nodes.
             */
            static std::size_t numPeriodicNodes(const Resolution& resolution) noexcept {
                return 3 * resolution[0] * resolution[1];
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
        };

    } // namespace grid

} // namespace monad
