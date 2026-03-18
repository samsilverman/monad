#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/element/hex8.hpp"

namespace monad {

    namespace grid {

        /**
         * @brief Topology rules for structured 3D Hex8 grids.
         */
        struct Hex8Topology {
            using Element = fem::Hex8;
            using Resolution = std::array<std::size_t, Element::Dim>;
            using Size = std::array<double, Element::Dim>;
            using Point = typename Element::Point;
            using ElementList = std::array<std::size_t, Element::NumNodes>;

            /**
             * @brief Number of nodes.
             *
             * @param[in] resolution Number of elements in each dimension.
             *
             * @returns Number of nodes.
             */
            static std::size_t numNodes(const Resolution& resolution) noexcept;

            /**
             * @brief Number of periodic nodes.
             *
             * @param[in] resolution Number of elements in each dimension.
             *
             * @returns Number of periodic nodes.
             */
            static std::size_t numPeriodicNodes(const Resolution& resolution) noexcept;

            /**
             * @brief Coordinates of a node.
             *
             * @param[in] index Node index.
             * @param[in] resolution Number of elements in each dimension.
             * @param[in] size Physical lengths in each dimension.
             *
             * @returns Coordinates of a node.
             */
            static Point node(std::size_t index, const Resolution& resolution, const Size& size) noexcept;

            /**
             * @brief Node indices of an element.
             *
             * @param[in] index Element index.
             * @param[in] resolution Number of elements in each dimension.
             *
             * @returns Node indices of an element.
             */
            static ElementList element(std::size_t index, const Resolution& resolution) noexcept;

            /**
             * @brief Periodic node indices of an element.
             *
             * @param[in] index Element index.
             * @param[in] resolution Number of elements in each dimension.
             *
             * @returns Periodic node indices of an element.
             */
            static ElementList periodicElement(std::size_t index, const Resolution& resolution) noexcept;
        };

    } // namespace grid

} // namespace monad
