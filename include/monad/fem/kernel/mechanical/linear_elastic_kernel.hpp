#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel_traits.hpp"
#include "monad/integration/integrate_matrix.hpp"
#include "monad/detail/eigen_utils.hpp"

namespace monad {

    namespace fem {

        namespace mechanical {

            /**
             * @brief Provides core FEM computations for a linear elastic element.
             *
             * This kernel implements the weak form of the linear elastic PDE:
             *
             * ∇·σ=∇·(Cε)=0
             *
             * where the displacements are decomposed into macroscopic and microscopic components:
             *
             * u=ū+ũ
             *
             * @tparam ElementT Element class (e.g. Quad4).
             */
            template <class Element>
            struct LinearElasticKernel {
                static_assert(Element::Dim == 2 || Element::Dim == 3, "Element spatial dimension must be 2 or 3.");

                /// @brief Number of dofs in the element.
                static constexpr int NumDofs = Element::Dim * Element::NumNodes;

                using Material = LinearElasticMaterial<Element::Dim>;
                using Traits = LinearElasticKernelTraits<Element::Dim>;

                using Point = typename Element::Point;
                using NodesMatrix = typename Element::NodesMatrix;
                using BMatrix = Eigen::Matrix<double, Material::VoigtSize, NumDofs>;

                /// @brief Element mechanical stiffness matrix type.
                using StiffnessMatrix = Eigen::Matrix<double, NumDofs, NumDofs>;

                /// @brief Element field matrix type.
                using FieldMatrix = Eigen::Matrix<double, NumDofs, Material::VoigtSize>;

                /**
                 * @brief Element B matrix evaluated at a local point.
                 *
                 * The B matrix relates an element's nodal displacements u to strains ε in Voigt
                 * notation at a local `point`:
                 *
                 * ε=Bu
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
                    const auto dNGlobal = J.inverse() * dN;

                    BMatrix B;
                    Traits::fillB(B, dNGlobal);

                    return B;
                }

                /**
                 * @brief Element mechanical stiffness matrix (left-hand side of the discretized weak form).
                 *
                 * Weak form lhs for an element e:
                 *
                 * Kₑ=∫_ΩₑBᵀCBdΩₑ
                 *
                 * @param[in] point Local point.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element mechanical stiffness matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
                 */
                static StiffnessMatrix lhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &C = material.materialTensor();

                    auto integrand = [&](const Point &point) -> StiffnessMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const BMatrix B = bMatrix(point, nodes);

                        return B.transpose() * C * B * J.determinant();
                    };

                    StiffnessMatrix K = integration::integrateMatrix(integrand, rule);

                    // Remove numerical artifacts
                    detail::symmetrize(K);

                    return K;
                }

                /**
                 * @brief Element force matrix (right-hand side of the discretized weak form).
                 *
                 * Weak form rhs for an element e:
                 *
                 * Fₑ=-∫_ΩₑBᵀCdΩₑε̄
                 *
                 * @param[in] point Local point.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element force matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
                 */
                static FieldMatrix rhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &C = material.materialTensor();

                    auto integrand = [&](const Point &point) -> FieldMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const BMatrix B = bMatrix(point, nodes);

                        return B.transpose() * C * J.determinant();
                    };

                    return -integration::integrateMatrix(integrand, rule);
                }
            };

        } // namespace mechanical

    } // namespace fem

} // namespace monad
