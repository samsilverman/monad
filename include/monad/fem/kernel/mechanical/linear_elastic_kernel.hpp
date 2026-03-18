#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/integration/integrate_matrix.hpp"
#include "monad/detail/eigen_utils.hpp"

namespace monad {

    namespace fem {

        namespace mechanical {

            /**
             * @brief Core FEM computations for a linear elastic element.
             *
             * This kernel implements the weak form of the linear elasticity:
             *
             * ```text
             * ∇·σ=∇·(Cε)=0
             * ```
             *
             * where the displacements u∈ℝᵈ are decomposed into macroscopic and
             * microscopic components:
             *
             * ```text
             * u=ū+ũ
             * ```
             *
             * @tparam ElementT Element type (e.g. Quad4).
             */
            template <class Element>
            struct LinearElasticKernel {
                static_assert(Element::Dim == 2 || Element::Dim == 3, "Element spatial dimension must be 2 or 3.");

                /// @brief Number of dofs in the element.
                static constexpr int NumDofs = Element::Dim * Element::NumNodes;

                using Material = material::LinearElasticMaterial<Element::Dim>;

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
                 * The B matrix maps nodal displacements u∈ℝᵈ to
                 * strains ε∈ℝᵛ at a local `point`:
                 *
                 * ```text
                 * ε=Bu
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
                    const auto dNGlobal = J.inverse() * dN;

                    BMatrix B = BMatrix::Zero();

                    if constexpr (Element::Dim == 2) {
                        // Normal strain ε₁₁
                        B(0, Eigen::seq(0, Eigen::indexing::last, 2)) = dNGlobal.row(0);

                        // Normal strain ε₂₂
                        B(1, Eigen::seq(1, Eigen::indexing::last, 2)) = dNGlobal.row(1);

                        // Shear strain ε₁₂
                        B(2, Eigen::seq(0, Eigen::indexing::last, 2)) = dNGlobal.row(1);
                        B(2, Eigen::seq(1, Eigen::indexing::last, 2)) = dNGlobal.row(0);
                    }
                    else {
                        // Normal strain ε₁₁
                        B(0, Eigen::seq(0, Eigen::indexing::last, 3)) = dNGlobal.row(0);

                        // Normal strain ε₂₂
                        B(1, Eigen::seq(1, Eigen::indexing::last, 3)) = dNGlobal.row(1);

                        // Normal strain ε₃₃
                        B(2, Eigen::seq(2, Eigen::indexing::last, 3)) = dNGlobal.row(2);

                        // Shear strain ε₁₂
                        B(3, Eigen::seq(0, Eigen::indexing::last, 3)) = dNGlobal.row(1);
                        B(3, Eigen::seq(1, Eigen::indexing::last, 3)) = dNGlobal.row(0);

                        // Shear strain ε₁₃
                        B(4, Eigen::seq(0, Eigen::indexing::last, 3)) = dNGlobal.row(2);
                        B(4, Eigen::seq(2, Eigen::indexing::last, 3)) = dNGlobal.row(0);

                        // Shear strain ε₂₃
                        B(5, Eigen::seq(1, Eigen::indexing::last, 3)) = dNGlobal.row(2);
                        B(5, Eigen::seq(2, Eigen::indexing::last, 3)) = dNGlobal.row(1);
                    }

                    return B;
                }

                /**
                 * @brief Element mechanical stiffness matrix evaluated at a local point.
                 *
                 * For an element e:
                 *
                 * ```text
                 * Kₑ=∫_ΩₑBᵀCBdΩₑ
                 * ```
                 *
                 * @param[in] material Linear elastic material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element mechanical stiffness matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
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

                    // Remove numerical asymmetry
                    detail::symmetrize(K);

                    return K;
                }

                /**
                 * @brief Element force matrix evaluated at a local point.
                 *
                 * For an element e:
                 *
                 * ```text
                 * Fₑ=-∫_ΩₑBᵀCdΩₑε̄
                 * ```
                 *
                 * @param[in] material Linear transport material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element force matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
                 */
                static FieldMatrix rhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &C = material.materialTensor();

                    auto integrand = [&](const Point &point) -> FieldMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const BMatrix B = bMatrix(point, nodes);

                        return B.transpose() * C * J.determinant();
                    };

                    // No need to multiply by ε̄=I for unit macroscopic strains
                    return -integration::integrateMatrix(integrand, rule);
                }
            };

        } // namespace mechanical

    } // namespace fem

} // namespace monad
