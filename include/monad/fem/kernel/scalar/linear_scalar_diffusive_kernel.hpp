#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/integration/integrate_matrix.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/material/transport/linear_transport_material.hpp"

namespace monad {

    namespace fem {

        namespace scalar {

            /// @brief Specifies the sign convention relating the physical field G to the scalar potential gradient ∇φ.
            enum class GradientConvention {
                /// @brief Physical field is defined as field G=−∇φ.
                Negative,
                /// @brief Physical field is defined as field G=∇φ.
                Positive
            };

            /**
             * @brief Provides core FEM computations for a linear scalar diffusive element.
             *
             * This kernel implements the weak form of the scalar diffusion PDE:
             *
             * ∇·J=∇·(K∇φ)=0
             *
             * where the scalar potentials are decomposed into macroscopic and microscopic components:
             *
             * φ=φ̄+φ̃
             *
             * @tparam ElementT Element class (e.g. Quad4).
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

                using Material = LinearTransportMaterial<Element::Dim>;

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
                 * The B matrix relates an element's nodal scalar potentials φ to scalar gradients at a
                 * local `point`:
                 *
                 * ∇φ=Bφ
                 *
                 * For a D-dimensional element:
                 *
                 * - D=2: point=[ξ η]ᵀ
                 *
                 * - D=3: point=[ξ η ζ]ᵀ
                 *
                 * @param[in] point Local point.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element B matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
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
                 * @brief Element diffusive stiffness matrix (left-hand side of the discretized weak form) evaluated at a local point.
                 *
                 * Weak form lhs for an element e:
                 *
                 * Kₑ=∫_ΩₑBᵀABdΩₑ
                 *
                 * @param[in] material Linear transport material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element diffusive stiffness matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
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

                    // Remove numerical artifacts
                    detail::symmetrize(K);

                    return K;
                }

                /**
                 * @brief Element source matrix (right-hand side of the discretized weak form) evaluated at a local point.
                 *
                 * Weak form rhs for an element e:
                 *
                 * Fₑ=-∫_ΩₑBᵀAdΩₑ∇φ̄
                 *
                 * @param[in] material Linear transport material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element source matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
                 */
                static FieldMatrix rhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &A = material.materialTensor();

                    auto integrand = [&](const Point &point) -> FieldMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const BMatrix B = bMatrix(point, nodes);

                        return B.transpose() * A * J.determinant();
                    };

                    return -integration::integrateMatrix(integrand, rule);
                }
            };

        } // namespace scalar

    } // namespace fem

} // namespace monad
