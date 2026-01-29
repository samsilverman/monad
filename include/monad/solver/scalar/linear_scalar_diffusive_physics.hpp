#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/kernel/scalar/linear_scalar_diffusive_kernel.hpp"
#include "monad/fem/operator/scalar/linear_scalar_diffusive_matrix_free_operator_traits.hpp"
#include "monad/solver/solver_options.hpp"

namespace monad {

    namespace solver {

        namespace scalar {

            /**
             * @brief Physics policy for linear scalar diffusive homogenization.
             *
             * @tparam Element Element class (e.g. Quad4).
             * @tparam C Gradient sign convention (GradientConvention::Negative or Positive).
             */
            template <class Element, fem::scalar::GradientConvention C>
            struct LinearScalarDiffusivePhysics {
                using Kernel = fem::scalar::LinearScalarDiffusiveKernel<Element, C>;
                using Material = typename Kernel::Material;

                /**
                 * @brief Number of macroscopic fields:
                 *
                 * - D=2: ∇φ̄₁₁, ∇φ̄₂₂.
                 *
                 * - D=3: ∇φ̄₁₁, ∇φ̄₂₂, ∇φ̄₃₃.
                 */
                static constexpr int NumMacroFields = Element::Dim;

                using OperatorTraits = fem::scalar::LinearScalarDiffusiveMatrixFreeOperatorTraits;

                /// @brief Global field matrix type.
                using FieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, NumMacroFields>;

                /// @brief Homogenized material tensor type.
                using MaterialTensor = typename Material::MaterialTensor;

                /**
                 * @brief Linear scalar diffusive homogenization results.
                 *
                 * @note Scalar potential fields are stored in arrays
                 * ordered by the prescribed macroscopic scalar field
                 * gradient loading directions:
                 *
                 * - D=2: ∇φ̄₁₁, ∇φ̄₂₂.
                 *
                 * - D=3: ∇φ̄₁₁, ∇φ̄₂₂, ∇φ̄₃₃.
                 *
                 * Each array entry corresponds to one independent macroscopic loading case.
                 */
                struct Results {
                    using NodalField = Eigen::VectorXd;

                    /// @brief Homogenized transport tensor K̄.
                    MaterialTensor KBar;

                    /// @brief Total nodal scalar potential fields φ=φ̄+φ̃.
                    std::array<NodalField, NumMacroFields> phi;

                    /// @brief Macroscopic nodal scalar potential fields φ̄.
                    std::array<NodalField, NumMacroFields> phiMacro;

                    /// @brief Microscopic nodal scalar potential fields φ̃.
                    std::array<NodalField, NumMacroFields> phiMicro;
                };

                /**
                 * @brief Macroscopic nodal scalar potentials.
                 *
                 * The macroscopic scalar potentials φ̄ are the scalar potentials
                 * resulting from macroscopic scalar potential gradients ∇φ̄:
                 *
                 * φ̄=(∇φ̄)·x
                 *
                 * @tparam Grid Grid class (e.g. Quad4Grid).
                 *
                 * @param[in] grid Periodic unit cell grid.
                 *
                 * @returns Macroscopic nodal scalar potentials.
                 */
                template <class Grid>
                static FieldMatrix macroscopicField(const GridBase<Grid, Element> &grid) {
                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = OperatorTraits::NumNodeDofs * numNodes;

                    FieldMatrix Phi = FieldMatrix::Zero(static_cast<int>(numDofs), NumMacroFields);

                    for (std::size_t i = 0; i < grid.numNodes(); ++i) {
                        Phi.row(static_cast<int>(i)) = Kernel::GradSign * grid.node(i);
                    }

                    return Phi;
                }

                /**
                 * @brief Constructs linear scalar diffusive homogenization results.
                 *
                 * @param[in] KBar Homogenized transport tensor.
                 * @param[in] Phi Total nodal scalar potential fields.
                 * @param[in] PhiMacro Macroscopic nodal scalar potential fields.
                 * @param[in] PhiMicro Microscopic nodal scalar potential fields.
                 * @param[in] fields Controls which nodal fields are stored in the solver results.
                 *
                 * @returns Linear scalar diffusive homogenization results.
                 */
                static Results makeResults(const MaterialTensor &KBar, const FieldMatrix &Phi, const FieldMatrix &PhiMacro, const FieldMatrix &PhiMicro, FieldSave fields) noexcept {
                    using NodalField = typename Results::NodalField;

                    Results results;

                    results.KBar = KBar;

                    if (wants(fields, FieldSave::Total)) {
                        std::array<NodalField, NumMacroFields> temp;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            temp[i] = Phi.col(static_cast<int>(i));
                        }

                        results.phi = temp;
                    }

                    if (wants(fields, FieldSave::Macro)) {
                        std::array<NodalField, NumMacroFields> temp;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            temp[i] = PhiMacro.col(static_cast<int>(i));
                        }

                        results.phiMacro = temp;
                    }

                    if (wants(fields, FieldSave::Micro)) {
                        std::array<NodalField, NumMacroFields> temp;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            temp[i] = PhiMicro.col(static_cast<int>(i));
                        }

                        results.phiMicro = temp;
                    }

                    return results;
                }
            };

        } // namespace scalar

    } // namespace solver

} // namespace monad
