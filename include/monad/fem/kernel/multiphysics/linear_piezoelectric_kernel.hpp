#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/fem/kernel/scalar/linear_scalar_diffusive_kernel.hpp"
#include "monad/integration/integrate_matrix.hpp"
#include "monad/detail/eigen_utils.hpp"

namespace monad {

    namespace fem {

        namespace multiphysics {

            /**
             * @brief Core FEM computations for a linear piezoelectric element.
             *
             * This kernel implements the weak form of the linear piezoelectricity:
             *
             * ```text
             * ∇·T=∇·(cS-dᵀE)=0
             * ∇·(D)=∇·(dS+ϵE)=0
             * ```
             *
             * where the displacements u∈ℝᵈ and electric potentials φ∈ℝ are decomposed into
             * macroscopic and microscopic components:
             *
             * ```text
             * u=ū+ũ
             * φ=φ̄+φ̃
             * ```
             *
             * @tparam Element Element type (e.g. Quad4).
             */
            template <class Element>
            struct LinearPiezoelectricKernel {
                static_assert(Element::Dim == 2 || Element::Dim == 3, "Element spatial dimension must be 2 or 3.");

                using MechanicalKernel = mechanical::LinearElasticKernel<Element>;
                using ElectricalKernel = scalar::LinearScalarDiffusiveKernel<Element, scalar::GradientConvention::Negative>;

                /// @brief Number of dofs in the element.
                static constexpr int NumDofs = ElectricalKernel::NumDofs + MechanicalKernel::NumDofs;

                using MechanicalMaterial = material::LinearElasticMaterial<Element::Dim>;
                using ElectricalMaterial = material::LinearTransportMaterial<Element::Dim>;
                using Material = material::LinearPiezoelectricMaterial<Element::Dim>;

                using Point = typename Element::Point;
                using NodesMatrix = typename Element::NodesMatrix;

                /// @brief Element coupling stiffness matrix type.
                using CouplingStiffnessMatrix = Eigen::Matrix<double, ElectricalKernel::NumDofs, MechanicalKernel::NumDofs>;

                /// @brief Element piezoelectric stiffness matrix type.
                using StiffnessMatrix = Eigen::Matrix<double, NumDofs, NumDofs>;

                /// @brief Element mechanical field matrix type induced by macroscopic electrical loading.
                using UPhiCouplingFieldMatrix = Eigen::Matrix<double, MechanicalKernel::NumDofs, Element::Dim>;

                /// @brief Element electrical field matrix type induced by macroscopic mechanical loading.
                using PhiUCouplingFieldMatrix = Eigen::Matrix<double, ElectricalKernel::NumDofs, Material::VoigtSize>;

                /// @brief Element electromechanical field matrix type.
                using FieldMatrix = Eigen::Matrix<double, NumDofs, Material::VoigtSize + Element::Dim>;

                /**
                 * @brief Element piezoelectric stiffness matrix evaluated at a local point.
                 *
                 * For an element e:
                 *
                 * ```text
                 * Kₑ = ⎡ (Kᵤᵤ)ₑ  -(Kᵤᵩ)ₑ⎤
                 *      ⎣-(Kᵩᵤ)ₑ  -(Kᵩᵩ)ₑ⎦
                 * ```
                 *
                 * - Element mechanical stiffness matrix:
                 *
                 * ```text
                 * (Kᵤᵤ)ₑ=∫_ΩₑBᵤᵀcBᵤdΩₑ
                 * ```
                 *
                 * - Element electrical stiffness matrix:
                 *
                 * ```text
                 * (Kᵩᵩ)ₑ=∫_ΩᵩBᵤᵀεBᵩdΩₑ
                 * ```
                 *
                 * - Element piezoelectric coupling stiffness matrix:
                 *
                 * ```text
                 * (Kᵩᵤ)ₑ=(Kᵤᵩ)ₑᵀ=∫_ΩᵩBᵩᵀdBᵤdΩₑ=(∫_ΩᵩBᵤᵀdᵀBᵩdΩₑ)ᵀ
                 * ```
                 *
                 * @param[in] material Linear piezoelectric material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element piezoelectric stiffness matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
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

                    const MechanicalMaterial mechanicalMaterial(material.stiffnessTensor());
                    const ElectricalMaterial electricalMaterial(material.permittivityTensor());

                    const auto Kuu = MechanicalKernel::lhs(mechanicalMaterial, nodes);
                    const auto Kphiphi = ElectricalKernel::lhs(electricalMaterial, nodes);
                    const CouplingStiffnessMatrix Kphiu = integration::integrateMatrix(integrand, rule);

                    StiffnessMatrix K;
                    K << Kuu, -Kphiu.transpose(),
                         -Kphiu, -Kphiphi;

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
                 * Fₑ = ⎡ (Fᵤ)ₑ⎤ = ⎡ (Fᵤᵤ)ₑ  (Fᵤᵩ)ₑ⎤
                 *      ⎣-(Fᵩ)ₑ⎦   ⎣-(Fᵩᵤ)ₑ -(Fᵩᵩ)ₑ⎦
                 * ```
                 *
                 * - (Fᵤ)ₑ=[-∫_ΩₑBᵤᵀcdΩₑT̄  ∫_ΩₑBᵤᵀdᵀdΩₑĒ]
                 *
                 * - (Fᵩ)ₑ=[-∫_ΩₑBᵩᵀdᵀdΩₑT̄ -∫_ΩₑBᵩᵀεdΩₑĒ]
                 *
                 * @param[in] material Linear piezoelectric material.
                 * @param[in] nodes Element nodes.
                 *
                 * @returns Element source matrix evaluated at `point`.
                 *
                 * @throws std::invalid_argument if `nodes` define a degenerate element.
                 * @throws std::invalid_argument if `nodes` define an inverted element.
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

                    const MechanicalMaterial mechanicalMaterial(material.stiffnessTensor());
                    const ElectricalMaterial electricalMaterial(material.permittivityTensor());

                    const auto Fuu = MechanicalKernel::rhs(mechanicalMaterial, nodes);
                    const auto Fphiphi = ElectricalKernel::rhs(electricalMaterial, nodes);

                    // No need to multiply by T̄=I for unit macroscopic strains
                    const PhiUCouplingFieldMatrix Fphiu = -integration::integrateMatrix(integrandPhiU, rule);

                    // No need to multiply by Ē=I for unit macroscopic electric fields
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
