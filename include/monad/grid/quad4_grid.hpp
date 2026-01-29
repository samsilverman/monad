#pragma once

#include "monad/grid/grid_2d_base.hpp"
#include "monad/fem/element/quad4.hpp"

namespace monad {

    /**
     * @brief 2D grid mesh of 4-node quadrilaterl (Quad4) elements.
     *
     * Example 3x2 grid:
     *
     * ```text
     *      o───────o───────o───────o
     *      │       │       │       │
     *      │       │       │       │
     *      │       │       │       │
     *      o───────o───────o───────o
     *      │       │       │       │
     *      │       │       │       │
     *      │       │       │       │
     *      o───────o───────o───────o
     * ```
     *
     * Global axes:
     *
     * - x runs left→right.
     *
     * - y runs bottom→top.
     */
    class Quad4Grid : public Grid2dBase<Quad4Grid, fem::Quad4> {
    public:
        using Grid2dBase<Quad4Grid, fem::Quad4>::Grid2dBase;

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
         *      8───────9───────10──────11
         *      │       │       │       │
         *      │       │       │       │
         *      │       │       │       │
         *      4───────5───────6───────7
         *      │       │       │       │
         *      │       │       │       │
         *      │       │       │       │
         *      0───────1───────2───────3
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
         *      8───────9───────10──────11
         *      │       │       │       │
         *      │   𝟹   │   𝟺   │   𝟻   │
         *      │       │       │       │
         *      4───────5───────6───────7
         *      │       │       │       │
         *      │   𝟶   │   𝟷   │   𝟸   │
         *      │       │       │       │
         *      0───────1───────2───────3
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
         *      0───────1───────2───────0
         *      │       │       │       │
         *      │   𝟹   │   𝟺   │   𝟻   │
         *      │       │       │       │
         *      3───────4───────5───────3
         *      │       │       │       │
         *      │   𝟶   │   𝟷   │   𝟸   │
         *      │       │       │       │
         *      0───────1───────2───────0
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
