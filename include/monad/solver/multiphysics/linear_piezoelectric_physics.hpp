#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/kernel/multiphysics/linear_piezoelectric_kernel.hpp"
#include "monad/fem/operator/multiphysics/linear_piezoelectric_matrix_free_operator_traits.hpp"
#include "monad/solver/mechanical/linear_elastic_physics.hpp"
#include "monad/solver/scalar/linear_scalar_diffusive_physics.hpp"
#include "monad/solver/solver_options.hpp"

namespace monad {

    namespace solver {

        namespace multiphysics {

            /**
             * @brief Physics policy for linear piezoelectric homogenization.
             *
             * @tparam Element Element class (e.g. Quad4).
             */
            template <class Element>
            struct LinearPiezoelectricPhysics {
                using Kernel = fem::multiphysics::LinearPiezoelectricKernel<Element>;
                using Material = typename Kernel::Material;

                /// @brief Spatial dimension (2 or 3).
                static constexpr int Dim = Element::Dim;

                /**
                 * @brief Number of macroscopic fields:
                 *
                 * - D=2: ε̄₁₁, ε̄₂₂, ε̄₁₂, Ē₁₁, Ē₂₂.
                 *
                 * - D=3: ε̄₁₁, ε̄₂₂, ε̄₃₃, ε̄₁₂, ε̄₁₃, ε̄₂₃, Ē₁₁, Ē₂₂, Ē₃₃.
                 */
                static constexpr int NumMacroFields = Material::VoigtSize + Dim;

                using OperatorTraits = fem::multiphysics::LinearPiezoelectricMatrixFreeOperatorTraits<Dim>;

                using MechanicalPhysics = solver::mechanical::LinearElasticPhysics<Element>;
                using ElectricalPhysics = solver::scalar::LinearScalarDiffusivePhysics<Element, fem::scalar::GradientConvention::Negative>;

                /// @brief Global electromechanical field matrix type.
                using FieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, NumMacroFields>;

                /// @brief Homogenized coupled constitutive operator type.
                using MaterialTensor = typename Material::MaterialTensor;

                /**
                 * @brief Results from linear piezoelectric homogenization.
                 *
                 * @note Displacement fields are stored in arrays
                 * ordered by the prescribed macroscopic strain
                 * loading directions:
                 *
                 * - D=2: ε̄₁₁, ε̄₂₂, ε̄₁₂, Ē₁₁, Ē₂₂.
                 *
                 * - D=3: ε̄₁₁, ε̄₂₂, ε̄₃₃, ε̄₁₂, ε̄₁₃, ε̄₂₃, Ē₁₁, Ē₂₂, Ē₃₃.
                 *
                 * Each array entry corresponds to one independent macroscopic loading case.
                 */
                struct Results {
                    using MechanicalNodalField = Eigen::Matrix<double, Eigen::Dynamic, Dim>;
                    using ElectricalNodalField = Eigen::VectorXd;

                    using LinearElasticMaterialTensor = typename MechanicalPhysics::Material::MaterialTensor;
                    using LinearDielectricMaterialTensor = typename ElectricalPhysics::Material::MaterialTensor;
                    using CouplingTensor = typename Material::CouplingTensor;

                    /// @brief Homogenized stiffness tensor c̄.
                    LinearElasticMaterialTensor cBar;

                    /// @brief Homogenized permittivity tensor ε̄.
                    LinearDielectricMaterialTensor epsilonBar;

                    /// @brief Homogenized piezoelectric coupling tensor d̄.
                    CouplingTensor dBar;

                    /// @brief Total nodal displacement fields u=ū+ũ.
                    std::array<MechanicalNodalField, NumMacroFields> u;

                    /// @brief Macroscopic nodal displacement fields ū.
                    std::array<MechanicalNodalField, NumMacroFields> uMacro;

                    /// @brief Microscopic nodal displacement fields ũ.
                    std::array<MechanicalNodalField, NumMacroFields> uMicro;

                    /// @brief Total nodal electric potential fields φ=φ̄+φ̃.
                    std::array<ElectricalNodalField, NumMacroFields> phi;

                    /// @brief Macroscopic nodal electric potential fields φ̄.
                    std::array<ElectricalNodalField, NumMacroFields> phiMacro;

                    /// @brief Microscopic nodal electric potential fields φ̃.
                    std::array<ElectricalNodalField, NumMacroFields> phiMicro;
                };

                /**
                 * @brief Macroscopic electromechanical fields.
                 *
                 * The macroscopic electromechanical fields X̄ are:
                 *
                 * ```text
                 * ⎡Ū 0⎤
                 * ⎣0 φ̄⎦
                 * ```
                 *
                 * - Ū are the displacements resulting from the macroscopic strains ε̄.
                 * - φ̄ are the electric potentials resulting from macroscopic electric fields Ē.
                 *
                 * @tparam Grid Grid class (e.g. Quad4Grid).
                 *
                 * @param[in] grid Periodic unit cell grid.
                 *
                 * @returns Macroscopic electromechanical fields.
                 */
                template <class Grid>
                static FieldMatrix macroscopicField(const GridBase<Grid, Element> &grid) {
                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = OperatorTraits::NumNodeDofs * numNodes;

                    FieldMatrix X = FieldMatrix::Zero(static_cast<int>(numDofs), NumMacroFields);

                    // Mechanical block
                    const auto U = MechanicalPhysics::macroscopicField(grid);
                    X.topLeftCorner(U.rows(), MechanicalPhysics::NumMacroFields) = U;

                    // Electrical block
                    const auto Phi = ElectricalPhysics::macroscopicField(grid);
                    X.bottomRightCorner(Phi.rows(), ElectricalPhysics::NumMacroFields) = Phi;

                    return X;
                }

                /**
                 * @brief Constructs linear piezoelectric homogenization results.
                 *
                 * @param[in] opBar Homogenized coupled constitutive operator.
                 * @param[in] X Total nodal fields.
                 * @param[in] XMacro Macroscopic nodal fields.
                 * @param[in] XMicro Microscopic nodal fields.
                 * @param[in] fields Controls which nodal fields are stored in the solver results.
                 *
                 * @returns Linear piezoelectric homogenization results.
                 */
                static Results makeResults(const MaterialTensor &opBar, const FieldMatrix &X, const FieldMatrix &XMacro, const FieldMatrix &XMicro, FieldSave fields) noexcept {
                    using MechanicalNodalField = typename Results::MechanicalNodalField;
                    using ElectricalNodalField = typename Results::ElectricalNodalField;

                    Results results;

                    results.cBar = opBar.topLeftCorner(Material::VoigtSize, Material::VoigtSize);
                    results.epsilonBar = -opBar.bottomRightCorner(Dim, Dim);
                    results.dBar = -opBar.bottomLeftCorner(Dim, Material::VoigtSize);

                    const auto numNodes = X.rows() / OperatorTraits::NumNodeDofs;
                    const auto numMechanicalDofs = numNodes * OperatorTraits::MechanicalTraits::NumNodeDofs;
                    const auto numElectricalDofs = numNodes * OperatorTraits::ElectricalTraits::NumNodeDofs;

                    if (wants(fields, FieldSave::Total)) {
                        const auto U = X.topRows(numMechanicalDofs);
                        const auto Phi = X.bottomRows(numElectricalDofs);

                        std::array<MechanicalNodalField, NumMacroFields> tempU;
                        std::array<ElectricalNodalField, NumMacroFields> tempPhi;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            tempU[i] = U.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                            tempPhi[i] = Phi.col(static_cast<int>(i));
                        }

                        results.u = tempU;
                        results.phi = tempPhi;
                    }

                    if (wants(fields, FieldSave::Macro)) {
                        const auto UMacro = XMacro.topRows(numMechanicalDofs);
                        const auto PhiMacro = XMacro.bottomRows(numElectricalDofs);

                        std::array<MechanicalNodalField, NumMacroFields> tempU;
                        std::array<ElectricalNodalField, NumMacroFields> tempPhi;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            tempU[i] = UMacro.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                            tempPhi[i] = PhiMacro.col(static_cast<int>(i));
                        }

                        results.uMacro = tempU;
                        results.phiMacro = tempPhi;
                    }

                    if (wants(fields, FieldSave::Micro)) {
                        const auto UMicro = XMicro.topRows(numMechanicalDofs);
                        const auto PhiMicro = XMicro.bottomRows(numElectricalDofs);

                        std::array<MechanicalNodalField, NumMacroFields> tempU;
                        std::array<ElectricalNodalField, NumMacroFields> tempPhi;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            tempU[i] = UMicro.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                            tempPhi[i] = PhiMicro.col(static_cast<int>(i));
                        }

                        results.uMicro = tempU;
                        results.phiMicro = tempPhi;
                    }

                    return results;
                }
            };

        } // namespace multiphysics

    } // namespace solver

} // namespace monad
