#pragma once

#include "monad/fem/element/element_base.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief 4-node quadrilateral (Quad4) finite element.
         *
         * The Quad4 element uses bilinear shape-functions on ξ,η∈[-1,1].
         *
         * Node ordering:
         *
         * ```text
         *      4───────3
         *      │       │
         *      │       │
         *      │       │
         *      1───────2
         * ```
         *
         * Local axes:
         *
         * - ξ runs left→right.
         *
         * - η runs bottom→top.
         */
        struct Quad4 : ElementBase<Quad4, 2, 4, 4> {
            /// @brief Local nodal coordinates.
            static NodesMatrix localNodes() noexcept;

            /**
             * @brief Shape functions evaluated at a local point.
             *
             * [ξ η]ᵀ → [N₁(ξ,η) ... N₄(ξ,η)]ᵀ
             *
             * @param[in] point Local point.
             *
             * @returns Shape functions evaluated at `point`.
             */
            static ShapeFuncVector shapeFunctions(const Point &point) noexcept;

            /**
             * @brief Shape function gradients evaluated at a local point.
             *
             * [ξ η]ᵀ →
             *
             * ```text
             * ⎡∂N₁/∂ξ ... ∂N₄/∂ξ⎤
             * ⎣∂N₁/∂η ... ∂N₄/∂η⎦
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
