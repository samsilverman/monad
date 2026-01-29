#pragma once

#include "monad/grid/grid_2d_base.hpp"
#include "monad/fem/element/quad8.hpp"

namespace monad {

    /**
     * @brief 2D grid mesh of 8-node quadrilaterl (Quad8) elements.
     *
     * Example 3x2 grid:
     *
     * ```text
     *      o───o───o───o───o───o───o
     *      │       │       │       │
     *      o       o       o       o
     *      │       │       │       │
     *      o───o───o───o───o───o───o
     *      │       │       │       │
     *      o       o       o       o
     *      │       │       │       │
     *      o───o───o───o───o───o───o
     * ```
     *
     * Global axes:
     *
     * - x runs left→right.
     *
     * - y runs bottom→top.
     */
    class Quad8Grid : public Grid2dBase<Quad8Grid, fem::Quad8> {
    public:
        using Grid2dBase<Quad8Grid, fem::Quad8>::Grid2dBase;

        /// @brief Number of nodes.
        std::size_t numNodes() const noexcept;

        /// @brief Number of periodic nodes.
        std::size_t numPeriodicNodes() const noexcept;

        /**
         * @brief Coordinates for a specific node.
         *
         * Example 3x2 grid:
         *
         * ```text
         *      8───18──9───19──10──20──11
         *      │       │       │       │
         *      25      26      27      28
         *      │       │       │       │
         *      4───15──5───16──6───17──7
         *      │       │       │       │
         *      21      22      23      24
         *      │       │       │       │
         *      0───12──1───13──2───14──3
         * ```
         *
         * @param[in] index Node index.
         *
         * @returns Coordinates.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numNodes()`).
         */
        Point node(std::size_t index) const;

        /**
         * @brief Node indices for a specific element.
         *
         * Example 3x2 grid:
         *
         * ```text
         *      8───18──9───19──10──20──11
         *      │       │       │       │
         *      25  𝟹   26  𝟺   27  𝟻   28
         *      │       │       │       │
         *      4───15──5───16──6───17──7
         *      │       │       │       │
         *      21  𝟶   22  𝟷   23  𝟸   24
         *      │       │       │       │
         *      0───12──1───13──2───14──3
         * ```
         *
         * @param[in] index Element index.
         *
         * @returns Node indices.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         */
        ElementList element(std::size_t index) const;

        /**
         * @brief Periodic node indices for a specific element.
         *
         * Example 3x2 grid:
         *
         * ```text
         *      0───6───1───7───2───8───0
         *      │       │       │       │
         *      15  𝟹   16  𝟺   17  𝟻   15
         *      │       │       │       │
         *      3───9───4───10──5───11──3
         *      │       │       │       │
         *      12  𝟶   13  𝟷   14  𝟸   13
         *      │       │       │       │
         *      0───6───1───7───2───8───0
         * ```
         *
         * @param[in] index Element index.
         *
         * @returns Periodic node indices.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         */
        ElementList periodicElement(std::size_t index) const;
    };

} // namespace monad
