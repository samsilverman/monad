#pragma once

#include "monad/fem/element/element_base.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief 20-node hexahedral (Hex20) finite element.
         *
         * The Hex20 element uses quadratic shape-functions on ξ,η,ζ∈[-1,1].
         *
         * Node ordering:
         *
         * ```text
         *           8────15───7
         *          /│        /│
         *        16 │      14 │
         *        / 20      /  19
         *       5────13───6   │
         *       │   │     │   │
         *       │   4───11│───3
         *      17  /      18 /
         *       │ 12      │ 10
         *       │/        │/
         *       1────9────2
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
        struct Hex20 : ElementBase<Hex20, 3, 20, 27> {
            /// @brief Local nodal coordinates.
            static NodesMatrix localNodes() noexcept;

            /**
             * @brief Shape functions evaluated at a local point.
             *
             * [ξ η ζ]ᵀ → [N₁(ξ,η,ζ) ... N₂₀(ξ,η,ζ)]ᵀ
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
             * ⎡∂N₁/∂ξ ... ∂N₂₀/∂ξ⎤
             * ⎪∂N₁/∂η ... ∂N₂₀/∂η⎪
             * ⎣∂N₁/∂ζ ... ∂N₂₀/∂ζ⎦
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
