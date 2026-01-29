#pragma once

#include "monad/fem/element/element_base.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief 8-node hexahedral (Hex8) finite element.
         *
         * The Hex8 element uses trilinear shape-functions on ξ,η,ζ∈[-1,1].
         *
         * Node ordering:
         *
         * ```text
         *           8─────────7
         *          /│        /│
         *         / │       / │
         *        /  │      /  │
         *       5─────────6   │
         *       │   │     │   │
         *       │   4─────│───3
         *       │  /      │  /
         *       │ /       │ /
         *       │/        │/
         *       1─────────2
         * ```
         *
         * Local axes:
         *
         * - ξ runs left→right.
         *
         * - η runs front→rear.
         *
         * - ζ runs bottom→top.
         */
        struct Hex8 : ElementBase<Hex8, 3, 8, 8> {
            /// @brief Local nodal coordinates.
            static NodesMatrix localNodes() noexcept;

            /**
             * @brief Shape functions evaluated at a local point.
             *
             * [ξ η ζ]ᵀ → [N₁(ξ,η,ζ) ... N₈(ξ,η,ζ)]ᵀ
             *
             * @param[in] point Local point.
             *
             * @returns Shape functions evaluated at `point`.
             */
            static ShapeFuncVector shapeFunctions(const Point &point) noexcept;

            /**
             * @brief Shape function gradients evaluated at a local point.
             *
             * [ξ η ζ]ᵀ →
             *
             * ```text
             * ⎡∂N₁/∂ξ ... ∂N₈/∂ξ⎤
             * ⎪∂N₁/∂η ... ∂N₈/∂η⎪
             * ⎣∂N₁/∂ζ ... ∂N₈/∂ζ⎦
             * ```
             *
             * @param[in] point Local point.
             *
             * @returns Shape function gradients evaluated at `point`.
             */
            static ShapeFuncGradMatrix gradShapeFunctions(const Point &point) noexcept;

            /// @brief Quadrature rule for integration.
            static QuadratureRule quadratureRule() noexcept;

            /// @brief Gmsh element type ID.
            static int gmshElementType() noexcept;

            /// @brief Gmsh element node ordering.
            static NodeIndicesList gmshNodeOrdering() noexcept;
        };

    } // namespace fem

} // namespace monad
