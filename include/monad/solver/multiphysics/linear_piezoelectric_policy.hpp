#pragma once

#include <array>
#include <cstddef>
#include <utility>
#include <Eigen/Core>
#include "monad/fem/kernel/multiphysics/linear_piezoelectric_kernel.hpp"
#include "monad/fem/operator/multiphysics/linear_piezoelectric_dof_traits.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/solver/mechanical/linear_elastic_policy.hpp"
#include "monad/solver/scalar/linear_scalar_diffusive_policy.hpp"
#include "monad/solver/solver_options.hpp"

namespace monad {

    namespace solver {

        namespace multiphysics {

            /**
             * @brief Physics policy for periodic linear piezoelectric homogenization.
             *
             * @tparam Element Element type (e.g. Quad4).
             */
            template <class Element>
            struct LinearPiezoelectricPolicy {
                using Kernel = fem::multiphysics::LinearPiezoelectricKernel<Element>;
                using DofTraits = fem::multiphysics::LinearPiezoelectricDofTraits<Element::Dim>;
                using Material = material::LinearPiezoelectricMaterial<Element::Dim>;
                using MaterialTensor = typename Material::MaterialTensor;
                using StiffnessTensor = typename Material::StiffnessTensor;
                using PermittivityTensor = typename Material::PermittivityTensor;
                using CouplingTensor = typename Material::CouplingTensor;
                using MechanicalPolicy = solver::mechanical::LinearElasticPolicy<Element>;
                using ElectricalPolicy = solver::scalar::LinearScalarDiffusivePolicy<Element, fem::scalar::GradientConvention::Negative>;

                /// @brief Spatial dimension (2 or 3).
                static constexpr int Dim = Element::Dim;

                /**
                 * @brief Number of macroscopic electromechanical load cases.
                 *
                 * The load cases are ordered as mechanical strain cases in
                 * Voigt order followed by electrical field cases in Cartesian
                 * component order:
                 *
                 * - 2D: S̄₁₁, S̄₂₂, S̄₁₂, Ē₁₁, Ē₂₂.
                 *
                 * - 3D: S̄₁₁, S̄₂₂, S̄₃₃, S̄₁₂, S̄₁₃, S̄₂₃, Ē₁₁, Ē₂₂, Ē₃₃.
                 */
                static constexpr int NumLoadCases = Material::VoigtSize + Dim;

                /// @brief Dof-major matrix type storing one global field per load case.
                using DofFieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, NumLoadCases>;

                /// @brief Nodal mechanical field type.
                using NodalFieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, Dim>;

                /// @brief Nodal electrical field type.
                using NodalFieldVector = Eigen::Vector<double, Eigen::Dynamic>;

                /// @brief List of nodal mechanical fields, one per load case.
                using MechanicalNodalFieldList = std::array<NodalFieldMatrix, NumLoadCases>;

                /// @brief List of nodal electrical fields, one per load case.
                using ElectricalNodalFieldList = std::array<NodalFieldVector, NumLoadCases>;

                /**
                 * @brief Linear piezoelectric homogenization results.
                 *
                 * The homogenized material tensors are always returned.
                 * Nodal fields are stored only when requested
                 * through `FieldSave`.
                 */
                struct Results {
                    /// @brief Homogenized coupled constitutive operator.
                    MaterialTensor opBar;

                    /// @brief Homogenized stiffness tensor.
                    StiffnessTensor cBar;

                    /// @brief Homogenized permittivity tensor.
                    PermittivityTensor epsilonBar;

                    /// @brief Homogenized piezoelectric coupling tensor.
                    CouplingTensor dBar;

                    /// @brief Total nodal displacements.
                    MechanicalNodalFieldList u;

                    /// @brief Macroscopic nodal displacements.
                    MechanicalNodalFieldList uMacro;

                    /// @brief Microscopic nodal displacements.
                    MechanicalNodalFieldList uMicro;

                    /// @brief Total nodal electric potentials.
                    ElectricalNodalFieldList phi;

                    /// @brief Macroscopic nodal electric potentials.
                    ElectricalNodalFieldList phiMacro;

                    /// @brief Microscopic nodal electric potentials.
                    ElectricalNodalFieldList phiMicro;
                };

                /**
                 * @brief Macroscopic electromechanical fields.
                 *
                 * The macroscopic field is assembled blockwise from:
                 *
                 * - Macroscopic displacements for unit macroscopic strains
                 *
                 * - Macroscopic electric potentials for unit macroscopic electric fields
                 *
                 * ```text
                 * ⎡Ū 0⎤
                 * ⎣0 Φ̄⎦
                 * ```
                 *
                 * @tparam Grid Grid type (e.g. Quad4Grid).
                 *
                 * @param[in] grid Grid.
                 *
                 * @returns Macroscopic electromechanical fields.
                 */
                template <class Grid>
                static DofFieldMatrix makeMacroscopicFields(const Grid &grid) noexcept {
                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = DofTraits::NumNodeDofs * numNodes;

                    DofFieldMatrix XBar = DofFieldMatrix::Zero(static_cast<int>(numDofs), NumLoadCases);

                    // Mechanical block
                    const auto UBar = MechanicalPolicy::makeMacroscopicFields(grid);
                    XBar.topLeftCorner(UBar.rows(), MechanicalPolicy::NumLoadCases) = UBar;

                    // Electrical block
                    const auto PhiBar = ElectricalPolicy::makeMacroscopicFields(grid);
                    XBar.bottomRightCorner(PhiBar.rows(), ElectricalPolicy::NumLoadCases) = PhiBar;

                    return XBar;
                }

                /**
                 * @brief Packages solver outputs into linear elastic results.
                 *
                 * @param[in] opBar Homogenized coupled constitutive operator.
                 * @param[in] X Total nodal fields.
                 * @param[in] XMacro Macroscopic nodal fields.
                 * @param[in] XMicro Microscopic nodal fields.
                 * @param[in] fields Controls which nodal fields are stored.
                 *
                 * @returns Homogenization results.
                 */
                static Results makeResults(const MaterialTensor &opBar, const DofFieldMatrix &X, const DofFieldMatrix &XMacro, const DofFieldMatrix &XMicro, FieldSave fields) noexcept {
                    Results results;
                    results.opBar = opBar;
                    results.cBar = opBar.topLeftCorner(Material::VoigtSize, Material::VoigtSize);
                    results.epsilonBar = -opBar.bottomRightCorner(Dim, Dim);
                    results.dBar = -opBar.bottomLeftCorner(Dim, Material::VoigtSize);

                    if (wants(fields, FieldSave::Total)) {
                        auto [u, phi] = makeStoredFields_(X);

                        results.u = u;
                        results.phi = phi;
                    }

                    if (wants(fields, FieldSave::Macro)) {
                        auto [uMacro, phiMacro] = makeStoredFields_(XMacro);

                        results.uMacro = uMacro;
                        results.phiMacro = phiMacro;
                    }

                    if (wants(fields, FieldSave::Micro)) {
                        auto [uMicro, phiMicro] = makeStoredFields_(XMicro);

                        results.uMicro = uMicro;
                        results.phiMicro = phiMicro;
                    }

                    return results;
                }

            private:
                /**
                 * @brief Splits a nodal coupled field block matrix into lists of nodal displacements and 
                 * and electric potentials.
                 *
                 * @param[in] X Nodal coupled field block matrix.
                 *
                 * @return Pair of lists for nodal displacements and nodal electric potentials.
                 */
                static std::pair<MechanicalNodalFieldList, ElectricalNodalFieldList> makeStoredFields_(const DofFieldMatrix &X) {
                    const int numNodes = static_cast<int>(X.rows()) / (Dim + 1);
                    const int numMechanicalDofs = numNodes * Dim;
                    const int numElectricalDofs = numNodes;

                    const auto U = X.topRows(numMechanicalDofs);
                    const auto Phi = X.bottomRows(numElectricalDofs);

                    MechanicalNodalFieldList outMechanical;
                    ElectricalNodalFieldList outElectrical;

                    for (std::size_t i = 0; i < NumLoadCases; ++i) {
                        const NodalFieldMatrix u = U.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                        const NodalFieldVector phi = Phi.col(static_cast<int>(i));

                        outMechanical[i] = u;
                        outElectrical[i] = phi;
                    }

                    return {outMechanical, outElectrical};
                }
            };

        } // namespace multiphysics

    } // namespace solver

} // namespace monad
