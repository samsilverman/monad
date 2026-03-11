#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/element/quad4.hpp"

namespace monad {

    namespace grid {

        /**
         * @brief Topology rules for 2D structured Quad4 grids.
         */
        struct Quad4Topology {
            using Element = fem::Quad4;
            using Resolution = std::array<std::size_t, Element::Dim>;
            using Size = std::array<double, Element::Dim>;
            using Point = typename Element::Point;
            using ElementList = std::array<std::size_t, Element::NumNodes>;

            /**
             * @brief Number of nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of nodes.
             */
            static std::size_t numNodes(const Resolution& resolution) noexcept {
                return (resolution[0] + 1) * (resolution[1] + 1);
            }

            /**
             * @brief Number of periodic nodes.
             *
             * @param[in] resolution How many cells in each dimension.
             *
             * @returns Number of periodic nodes.
             */
            static std::size_t numPeriodicNodes(const Resolution& resolution) noexcept {
                return resolution[0] * resolution[1];
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
                const std::size_t i = index % (nx + 1);
                const std::size_t j = index / (nx + 1);

                const double dx = size[0] / static_cast<double>(nx);
                const double dy = size[1] / static_cast<double>(resolution[1]);

                const double x = static_cast<double>(i) * dx;
                const double y = static_cast<double>(j) * dy;

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
        };

    } // namespace grid

} // namespace monad
