#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/fem/kernel/scalar/linear_scalar_diffusive_kernel.hpp"
#include "monad/integration/integrate_matrix.hpp"

namespace monad {

    namespace fem {

        namespace multiphysics {

            /**
             * @brief Provides core FEM computations for a linear piezoelectric element.
             *
             * This kernel implements the weak form of the linear piezoelectric PDE:
             *
             * ∇·S=∇·(cT-dᵀE)=0
             * ∇·(-D)=∇·(-dT-εE)=0
             *
             * where the displacements and electric potentials are decomposed into macroscopic and microscopic components:
             *
             * u=ū+ũ
             * φ=φ̄+φ̃
             *
             * @tparam Element Element class (e.g. Quad4).
             */
            template <class Element>
            struct LinearPiezoelectricKernel {
                static_assert(Element::Dim == 2 || Element::Dim == 3, "Element spatial dimension must be 2 or 3.");

                using MechanicalKernel = mechanical::LinearElasticKernel<Element>;
                using ElectricalKernel = scalar::LinearScalarDiffusiveKernel<Element, scalar::GradientConvention::Negative>;

                /// @brief Number of dofs in the element.
                static constexpr int NumDofs = ElectricalKernel::NumDofs + MechanicalKernel::NumDofs;

                using Material = LinearPiezoelectricMaterial<LinearElasticMaterial<Element::Dim>, LinearTransportMaterial<Element::Dim>>;

                using Point = typename Element::Point;
                using NodesMatrix = typename Element::NodesMatrix;

                /// @brief Element coupling stiffness matrix type.
                using CouplingStiffnessMatrix = Eigen::Matrix<double, ElectricalKernel::NumDofs, MechanicalKernel::NumDofs>;

                /// @brief Element piezoelectric stiffness matrix type.
                using StiffnessMatrix = Eigen::Matrix<double, NumDofs, NumDofs>;

                /// @brief Element mechanical field (induced by macroscopic electrical loading) matrix type.
                using UPhiCouplingFieldMatrix = Eigen::Matrix<double, MechanicalKernel::NumDofs, Element::Dim>;

                /// @brief Element electrical field (induced by macroscopic mechanical loading) matrix type.
                using PhiUCouplingFieldMatrix = Eigen::Matrix<double, ElectricalKernel::NumDofs, Material::VoigtSize>;

                /// @brief Element electromechanical field matrix type.
                using FieldMatrix = Eigen::Matrix<double, NumDofs, Material::VoigtSize + Element::Dim>;

                /**
                 * @brief Element piezoelectric stiffness matrix (left-hand side of the discretized weak form).
                 *
                 * Weak form lhs for an element e:
                 *
                 * ```text
                 * Kₑ = ⎡ (Kᵤᵤ)ₑ  -(Kᵤᵩ)ₑ⎤
                 *      ⎣-(Kᵩᵤ)ₑ  -(Kᵩᵩ)ₑ⎦
                 * ```
                 *
                 * - (Kᵤᵤ)ₑ=∫_ΩₑBᵤᵀsBᵤdΩₑ (element mechanical stiffness matrix)
                 *
                 * - (Kᵩᵩ)ₑ=∫_ΩᵩBᵤᵀεBᵩdΩₑ (element electrical stiffness matrix)
                 *
                 * - (Kᵩᵤ)ₑ=(Kᵤᵩ)ₑᵀ=∫_ΩᵩBᵩᵀdBᵤdΩₑ=(∫_ΩᵩBᵤᵀdᵀBᵩdΩₑ)ᵀ (element coupling stiffness matrix)
                 *
                 * @param[in] point Local point.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element piezoelectric stiffness matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
                 */
                static StiffnessMatrix lhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &d = material.couplingTensor();

                    auto integrand = [&](const Point &point) -> CouplingStiffnessMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const auto Bu = MechanicalKernel::bMatrix(point, nodes);
                        const auto Bphi = ElectricalKernel::bMatrix(point, nodes);

                        return Bphi.transpose() * d * Bu * J.determinant();
                    };

                    const auto Kuu = MechanicalKernel::lhs(material.elasticMaterial(), nodes);
                    const auto Kphiphi = ElectricalKernel::lhs(material.dielectricMaterial(), nodes);
                    const CouplingStiffnessMatrix Kphiu = integration::integrateMatrix(integrand, rule);

                    StiffnessMatrix K;
                    K << Kuu, -Kphiu.transpose(),
                         -Kphiu, -Kphiphi;

                    return K;
                }

                /**
                 * @brief Element source matrix (right-hand side of the discretized weak form).
                 *
                 * Weak form rhs for an element e:
                 *
                 * ```text
                 * Fₑ = ⎡ (Fᵤ)ₑ⎤ = ⎡ (Fᵤᵤ)ₑ  (Fᵤᵩ)ₑ⎤
                 *      ⎣-(Fᵩ)ₑ⎦   ⎣-(Fᵩᵤ)ₑ -(Fᵩᵩ)ₑ⎦
                 * ```
                 *
                 * - (Fᵤ)ₑ=[-∫_ΩₑBᵤᵀsdΩₑT̄  ∫_ΩₑBᵤᵀdᵀdΩₑĒ]
                 *
                 * - (Fᵩ)ₑ=[-∫_ΩₑBᵩᵀdᵀdΩₑT̄ -∫_ΩₑBᵩᵀεdΩₑĒ]
                 *
                 * @param[in] point Local point.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element source matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate or inverted element geometry.
                 */
                static FieldMatrix rhs(const Material &material, const NodesMatrix &nodes) {
                    const auto rule = Element::quadratureRule();

                    const auto &d = material.couplingTensor();

                    auto integrandPhiU = [&](const Point &point) -> PhiUCouplingFieldMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const auto B = ElectricalKernel::bMatrix(point, nodes);

                        return B.transpose() * d * J.determinant();
                    };

                    auto integrandUPhi = [&](const Point &point) -> UPhiCouplingFieldMatrix {
                        const auto J = Element::jacobian(point, nodes);
                        const auto B = MechanicalKernel::bMatrix(point, nodes);

                        return B.transpose() * d.transpose() * J.determinant();
                    };

                    const auto Fuu = MechanicalKernel::rhs(material.elasticMaterial(), nodes);
                    const auto Fphiphi = ElectricalKernel::rhs(material.dielectricMaterial(), nodes);
                    const PhiUCouplingFieldMatrix Fphiu = -integration::integrateMatrix(integrandPhiU, rule);
                    const UPhiCouplingFieldMatrix Fuphi = integration::integrateMatrix(integrandUPhi, rule);

                    FieldMatrix F;
                    F << Fuu, Fuphi,
                         -Fphiu, -Fphiphi;

                    return F;
                }
            };

        } // namespace multiphysics

    } // namespace fem

} // namespace monad
