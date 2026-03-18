#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/integration/integrate_matrix.hpp"
#include "monad/detail/eigen_utils.hpp"

namespace monad {

    namespace fem {

        namespace scalar {

            /// @brief Sign convention for the derived field ∇φ.
            enum class GradientConvention {
                /// @brief Derived field is defined as −∇φ.
                Negative,
                /// @brief Derived field is defined as ∇φ.
                Positive
            };

            /**
             * @brief Core FEM computations for a linear scalar diffusive element.
             *
             * This kernel implements the weak form of scalar diffusion:
             *
             * ```text
             * ∇·J=∇·(K∇φ)=0
             * ```
             *
             * where the scalar potentials φ∈ℝ are decomposed into macroscopic and
             * microscopic components:
             *
             * ```text
             * φ=φ̄+φ̃
             * ```
             *
             * @tparam ElementT Element type (e.g. Quad4).
             * @tparam C Gradient sign convention (GradientConvention::Negative or Positive).
             */
            template <class ElementT, GradientConvention C>
            struct LinearScalarDiffusiveKernel {
                static_assert(ElementT::Dim == 2 || ElementT::Dim == 3, "Element spatial dimension must be 2 or 3.");

                using Element = ElementT;

                /// @brief Number of dofs in the element.
                static constexpr int NumDofs = Element::NumNodes;

                /// @brief Gradient sign.
                static constexpr double GradSign = (C == GradientConvention::Negative ? -1.0 : 1.0);

                using Material = material::LinearTransportMaterial<Element::Dim>;

                using Point = typename Element::Point;
                using NodesMatrix = typename Element::NodesMatrix;
                using BMatrix = Eigen::Matrix<double, Element::Dim, Element::NumNodes>;

                /// @brief Element diffusive stiffness matrix type.
                using StiffnessMatrix = Eigen::Matrix<double, Element::NumNodes, Element::NumNodes>;

                /// @brief Element field matrix type.
                using FieldMatrix = Eigen::Matrix<double, Element::NumNodes, Element::Dim>;

                /**
                 * @brief Element B matrix evaluated at a local point.
                 *
                 * The B matrix masp nodal scalar potentials φ∈ℝ to
                 * scalar potential gradients ∇φ∈ℝᵈ at a local `point`:
                 *
                 * ```text
                 * ∇φ=Bφ
                 * ```
                 *
                 * @param[in] point Local point.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element B matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
                 */
                static BMatrix bMatrix(const Point &point, const NodesMatrix &nodes) {
                    const auto J = Element::jacobian(point, nodes);
                    const double detJ = J.determinant();

                    if (detJ == 0) {
                        throw std::invalid_argument("Element Jacobian determinant is zero. The element is degenerate.");
                    }

                    if (detJ < 0) {
                        throw std::invalid_argument("Element Jacobian determinant is negative. The element is inverted.");
                    }

                    const auto dN = Element::gradShapeFunctions(point);

                    return GradSign * J.inverse() * dN;
                }

                /**
                 * @brief Element diffusive stiffness matrix evaluated at a local point.
                 *
                 * For an element e:
                 *
                 * ```text
                 * Kₑ=∫_ΩₑBᵀKBdΩₑ
                 * ```
                 *
                 * @param[in] material Linear transport material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element diffusive stiffness matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
                 */
                static StiffnessMatrix lhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &A = material.materialTensor();

                    auto integrand = [&](const Point &point) -> StiffnessMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const BMatrix B = bMatrix(point, nodes);

                        return B.transpose() * A * B * J.determinant();
                    };

                    StiffnessMatrix K = integration::integrateMatrix(integrand, rule);

                    // Remove numerical asymmetry
                    detail::symmetrize(K);

                    return K;
                }

                /**
                 * @brief Element source matrix evaluated at a local point.
                 *
                 * For an element e:
                 *
                 * ```text
                 * Fₑ=-∫_ΩₑBᵀAdΩₑ∇φ̄
                 * ```
                 *
                 * @param[in] material Linear transport material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element source matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
                 */
                static FieldMatrix rhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &A = material.materialTensor();

                    auto integrand = [&](const Point &point) -> FieldMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const BMatrix B = bMatrix(point, nodes);

                        return B.transpose() * A * J.determinant();
                    };

                    // No need to multiply by ∇φ̄=I for unit macroscopic scalar potential gradients
                    return -integration::integrateMatrix(integrand, rule);
                }
            };

        } // namespace scalar

    } // namespace fem

} // namespace monad
