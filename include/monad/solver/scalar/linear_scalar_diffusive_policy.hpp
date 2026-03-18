#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/kernel/scalar/linear_scalar_diffusive_kernel.hpp"
#include "monad/fem/operator/scalar/linear_scalar_diffusive_dof_traits.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/solver/solver_options.hpp"

namespace monad {

    namespace solver {

        namespace scalar {

            /**
             * @brief Physics policy for periodic linear scalar diffusive homogenization.
             *
             * @tparam Element Element type (e.g. Quad4).
             * @tparam C Gradient sign convention (GradientConvention::Negative or Positive).
             */
            template <class Element, fem::scalar::GradientConvention C>
            struct LinearScalarDiffusivePolicy {
                using Kernel = fem::scalar::LinearScalarDiffusiveKernel<Element, C>;
                using DofTraits = fem::scalar::LinearScalarDiffusiveDofTraits;
                using Material = material::LinearTransportMaterial<Element::Dim>;
                using MaterialTensor = typename Material::MaterialTensor;

                /// @brief Spatial dimension (2 or 3).
                static constexpr int Dim = Element::Dim;

                /**
                 * @brief Number of macroscopic scalar potential gradient load cases.
                 *
                 * The load cases are ordered in Cartesian component order:
                 *
                 * - 2D: ∇φ̄₁₁, ∇φ̄₂₂.
                 *
                 * - 3D: ∇φ̄₁₁, ∇φ̄₂₂, ∇φ̄₃₃.
                 */
                static constexpr int NumLoadCases = Dim;

                /// @brief Dof-major matrix type storing one global field per load case.
                using DofFieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, NumLoadCases>;

                /// @brief Nodal field type.
                using NodalFieldVector = Eigen::Vector<double, Eigen::Dynamic>;

                /// @brief List of nodal fields, one per load case.
                using NodalFieldList = std::array<NodalFieldVector, NumLoadCases>;

                /**
                 * @brief Linear scalar diffusive homogenization results.
                 *
                 * The homogenized transport tensor is always returned.
                 * Nodal scalar potentials are stored only when requested
                 * through `FieldSave`.
                 *
                 * Each array entry corresponds to one independent macroscopic loading case.
                 */
                struct Results {
                    /// @brief Homogenized transport tensor.
                    MaterialTensor KBar;

                    /// @brief Total nodal scalar potential fields.
                    NodalFieldList phi;

                    /// @brief Macroscopic nodal scalar potential fields.
                    NodalFieldList phiMacro;

                    /// @brief Microscopic nodal scalar potential fields.
                    NodalFieldList phiMicro;
                };

                /**
                 * @brief Macroscopic nodal scalar potentials.
                 *
                 * For each prescribed macroscopic scalar potential gradient ∇φ̄,
                 * the macroscopic scalar potentials are:
                 *
                 * ```text
                 * φ̄=(∇φ̄)·x
                 * ```
                 *
                 * @tparam Grid Grid type (e.g. Quad4Grid).
                 *
                 * @param[in] grid Grid.
                 *
                 * @returns Macroscopic nodal scalar potentials.
                 */
                template <class Grid>
                static DofFieldMatrix makeMacroscopicFields(const Grid &grid) {
                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = DofTraits::NumNodeDofs * numNodes;

                    DofFieldMatrix PhiBar = DofFieldMatrix::Zero(static_cast<int>(numDofs), NumLoadCases);

                    for (std::size_t i = 0; i < numNodes; ++i) {
                        PhiBar.row(static_cast<int>(i)) = Kernel::GradSign * grid.node(i);
                    }

                    return PhiBar;
                }

                /**
                 * @brief Packages solver outputs into linear elastic results.
                 *
                 * @param[in] KBar Homogenized transport tensor.
                 * @param[in] Phi Total nodal scalar potentials.
                 * @param[in] PhiMacro Macroscopic nodal scalar potentials.
                 * @param[in] PhiMicro Microscopic nodal scalar potentials.
                 * @param[in] fields Controls which nodal fields are stored.
                 *
                 * @returns Homogenization results.
                 */
                static Results makeResults(const MaterialTensor &KBar, const DofFieldMatrix &Phi, const DofFieldMatrix &PhiMacro, const DofFieldMatrix &PhiMicro, FieldSave fields) noexcept {
                    Results results;
                    results.KBar = KBar;

                    if (wants(fields, FieldSave::Total)) {
                        results.phi = makeStoredFields_(Phi);
                    }

                    if (wants(fields, FieldSave::Macro)) {
                        results.phiMacro = makeStoredFields_(PhiMacro);
                    }

                    if (wants(fields, FieldSave::Micro)) {
                        results.phiMicro = makeStoredFields_(PhiMicro);
                    }

                    return results;
                }
            
            private:
                /**
                 * @brief Converts a nodal scalar potential block matrix into a list of nodal scalar potentials.
                 *
                 * @param[in] Phi Nodal scalar potential block matrix.
                 *
                 * @return List of nodal scalar potentials.
                 */
                static NodalFieldList makeStoredFields_(const DofFieldMatrix &Phi) {
                    NodalFieldList out;

                    for (std::size_t i = 0; i < NumLoadCases; ++i) {
                        const NodalFieldVector phi = Phi.col(static_cast<int>(i));
                        out[i] = phi;
                    }

                    return out;
                }
            };

        } // namespace scalar

    } // namespace solver

} // namespace monad
