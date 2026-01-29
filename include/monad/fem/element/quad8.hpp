#pragma once

#include "monad/fem/element/element_base.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief 8-node quadrilateral (Quad8) finite element.
         *
         * The Quad8 element uses quadratic shape-functions on ξ,η∈[-1,1].
         *
         * Node ordering:
         *
         * ```text
         *      4───7───3
         *      │       │
         *      8       6
         *      │       │
         *      1───5───2
         * ```
         *
         * Local axes:
         *
         * - ξ runs left→right.
         *
         * - η runs bottom→top.
         */
        struct Quad8 : ElementBase<Quad8, 2, 8, 9> {
            /// @brief Local nodal coordinates.
            static NodesMatrix localNodes() noexcept;

            /**
             * @brief Shape functions evaluated at a local point.
             *
             * [ξ η]ᵀ → [N₁(ξ,η) ... N₈(ξ,η)]ᵀ
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
             * ⎡∂N₁/∂ξ ... ∂N₈/∂ξ⎤
             * ⎣∂N₁/∂η ... ∂N₈/∂η⎦
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
