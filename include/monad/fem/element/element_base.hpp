#pragma once

#include <array>
#include <cstddef>
#include <cmath>
#include <Eigen/Dense>
#include "monad/integration/quadrature_rule.hpp"
#include "monad/integration/integrate_scalar.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief Curiously recurring template pattern (CRTP) base class for FEM elements.
         *
         * @tparam Derived Concrete element class (for CRTP).
         * @tparam D Spatial dimension (2 or 3).
         * @tparam K Number of nodes in the element.
         * @tparam N Number of integration points in the element.
         *
         * @note `Derived` classes must provide:
         *
         * - `static NodesMatrix localNodes() noexcept`
         *
         * - `static ShapeFuncVector shapeFunctions(const Point &point) noexcept`
         *
         * - `static ShapeFuncGradMatrix gradShapeFunctions(const Point &point) noexcept`
         *
         * - `static QuadratureRule quadratureRule() noexcept`
         *
         * - `static int gmshElementType() noexcept`
         *
         * - `static NodeIndicesList gmshNodeOrdering() noexcept`
         */
        template <class Derived, int D, int K, int N>
        struct ElementBase {
            static_assert(D == 2 || D == 3, "Spatial dimension must be 2 or 3.");
            static_assert(K > 0, "Number of nodes must be positive.");
            static_assert(N > 0, "Number of integration points must be positive.");

            /// @brief Spatial dimension (2 or 3).
            static constexpr int Dim = D;

            /// @brief Number of nodes.
            static constexpr int NumNodes = K;

            /// @brief Number of integration points.
            static constexpr int NumIntegrationPoints = N;

            using Point = Eigen::Vector<double, Dim>;
            using NodesMatrix = Eigen::Matrix<double, NumNodes, Dim>;
            using ShapeFuncVector = Eigen::Vector<double, NumNodes>;
            using ShapeFuncGradMatrix = Eigen::Matrix<double, Dim, NumNodes>;
            using JacobianMatrix = Eigen::Matrix<double, Dim, Dim>;

            using QuadratureRule = integration::QuadratureRule<Dim, NumIntegrationPoints>;

            using NodeIndicesList = std::array<std::size_t, NumNodes>;

            /// @brief Local nodal coordinates.
            static NodesMatrix localNodes() noexcept {
                return Derived::localNodes();
            }

            /**
             * @brief Shape functions evaluated at a local point.
             *
             * For a D-dimensional element:
             *
             * - D=2: [ξ η]ᵀ → [N₁(ξ,η) ... Nₖ(ξ,η)]ᵀ
             *
             * - D=3: [ξ η ζ]ᵀ → [N₁(ξ,η,ζ) ... Nₖ(ξ,η,ζ)]ᵀ
             *
             * @param[in] point Local point.
             *
             * @returns Shape functions evaluated at `point`.
             */
            static ShapeFuncVector shapeFunctions(const Point &point) noexcept {
                return Derived::shapeFunctions(point);
            }

            /**
             * @brief Shape function gradients evaluated at a local point.
             *
             * For a D-dimensional element:
             *
             * - D=2: [ξ η]ᵀ →
             *
             * ```text
             * ⎡∂N₁/∂ξ ... ∂Nₖ/∂ξ⎤
             * ⎣∂N₁/∂η ... ∂Nₖ/∂η⎦
             * ```
             *
             * - D=3: [ξ η ζ]ᵀ →
             *
             * ```text
             * ⎡∂N₁/∂ξ ... ∂Nₖ/∂ξ⎤
             * ⎪∂N₁/∂η ... ∂Nₖ/∂η⎪
             * ⎣∂N₁/∂ζ ... ∂Nₖ/∂ζ⎦
             * ```
             *
             * @param[in] point Local point.
             *
             * @returns Shape function gradients evaluated at `point`.
             */
            static ShapeFuncGradMatrix gradShapeFunctions(const Point &point) noexcept {
                return Derived::gradShapeFunctions(point);
            }

            /**
             * @brief Jacobian matrix evaluated at a local point.
             *
             * For a D-dimensional element:
             *
             * - D=2: [ξ η]ᵀ →
             *
             * ```text
             * ⎡∂x/∂ξ ∂x/∂η⎤
             * ⎣∂y/∂ξ ∂y/∂η⎦
             * ```
             *
             * - D=3: [ξ η ζ]ᵀ →
             *
             * ```text
             * ⎡∂x/∂ξ ∂x/∂η ∂x/∂ζ⎤
             * ⎪∂y/∂ξ ∂y/∂η ∂y/∂ζ⎪
             * ⎣∂z/∂ξ ∂z/∂η ∂z/∂ζ⎦
             * ```
             *
             * @param[in] point Local point.
             * @param[in] nodes Element nodes.
             *
             * @returns Jacobian matrix evaluated at `point`.
             */
            static JacobianMatrix jacobian(const Point &point, const NodesMatrix &nodes) noexcept {
                return gradShapeFunctions(point) * nodes;
            }

            /// @brief Quadrature rule for integration.
            static QuadratureRule quadratureRule() noexcept {
                return Derived::quadratureRule();
            }

            /**
             * @brief Area (2D) or volume (3D).
             *
             * Computes area/volume via:
             *
             * ∫_Ω|det(J)|dΩ
             *
             * @param[in] nodes Element nodes.
             *
             * @returns Area (2D) or volume (3D).
             */
            static double measure(const NodesMatrix &nodes) noexcept {
                auto integrand = [&](const Point &point) -> double {
                    const JacobianMatrix J = jacobian(point, nodes);
                    return std::abs(J.determinant());
                };

                return integration::integrateScalar(integrand, quadratureRule());
            }

            /// @brief Gmsh element type ID.
            static int gmshElementType() noexcept {
                return Derived::gmshElementType();
            }

            /// @brief Gmsh element node ordering.
            static NodeIndicesList gmshNodeOrdering() noexcept {
                return Derived::gmshNodeOrdering();
            }
        };

    } // namespace fem

} // namespace monad
