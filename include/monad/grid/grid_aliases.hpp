#pragma once

#include "monad/grid/structured_grid.hpp"
#include "monad/grid/quad4_topology.hpp"
#include "monad/grid/quad8_topology.hpp"
#include "monad/grid/hex8_topology.hpp"
#include "monad/grid/hex20_topology.hpp"

namespace monad {

    /**
     * @brief 2D structured grid of 4-node quadrilateral (Quad4) elements.
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
    using Quad4Grid = grid::StructuredGrid<grid::Quad4Topology>;

    /**
     * @brief 2D structured grid of 8-node quadrilateral (Quad8) elements.
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
    using Quad8Grid = grid::StructuredGrid<grid::Quad8Topology>;

        /**
     * @brief 3D structured grid of 8-node hexahedral (Hex8) elements.
     *
     * Example 3x2x2 grid:
     *
     * ```text
     *               o─────────o─────────o─────────o
     *              /         /         /         /│
     *             /         /         /         / │
     *            /         /         /         /  │
     *           o─────────o─────────o─────────o   │
     *          /         /         /         /│   │
     *         /         /         /         / │   o
     *        /         /         /         /  │  /│
     *       o─────────o─────────o─────────o   │ / │
     *       │         │         │         │   │/  │
     *       │         │         │         │   o   │
     *       │         │         │         │  /│   │
     *       │         │         │         │ / │   o
     *       │         │         │         │/  │  /
     *       o─────────o─────────o─────────o   │ /
     *       │         │         │         │   │/
     *       │         │         │         │   o
     *       │         │         │         │  /
     *       │         │         │         │ /
     *       │         │         │         │/
     *       o─────────o─────────o─────────o
     * ```
     *
     * Global axes:
     *
     * - x runs left→right.
     *
     * - y runs front→back.
     *
     * - z runs bottom→top.
     */
    using Hex8Grid = grid::StructuredGrid<grid::Hex8Topology>;

    /**
     * @brief 3D structured grid of 20-node hexahedral (Hex20) elements.
     *
     * Example 3x2x2 grid:
     *
     * ```text
     *               o────o────o────o────o────o────o
     *              /         /         /         /│
     *             o         o         o         o │
     *            /         /         /         /  o
     *           o────o────o────o────o────o────o   │
     *          /         /         /         /│   │
     *         o         o         o         o │   o
     *        /         /         /         /  o  /│
     *       o────o────o────o────o────o────o   │ o │
     *       │         │         │         │   │/  o
     *       │         │         │         │   o   │
     *       o         o         o         o  /│   │
     *       │         │         │         │ o │   o
     *       │         │         │         │/  o  /
     *       o────o────o────o────o────o────o   │ o
     *       │         │         │         │   │/
     *       │         │         │         │   o
     *       o         o         o         o  /
     *       │         │         │         │ o
     *       │         │         │         │/
     *       o────o────o────o────o────o────o
     * ```
     *
     * Global axes:
     *
     * - x runs left→right.
     *
     * - y runs front→back.
     *
     * - z runs bottom→top.
     */
    using Hex20Grid = grid::StructuredGrid<grid::Hex20Topology>;

} // namespace monad
