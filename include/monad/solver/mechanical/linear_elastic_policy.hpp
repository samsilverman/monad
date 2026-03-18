#pragma once

#include <array>
#include <cstddef>
#include <Eigen/Core>
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/fem/operator/mechanical/linear_elastic_dof_traits.hpp"
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/solver/solver_options.hpp"

namespace monad {

    namespace solver {

        namespace mechanical {

            /**
             * @brief Physics policy for periodic linear elastic homogenization.
             *
             * @tparam Element Element type (e.g. Quad4).
             */
            template <class Element>
            struct LinearElasticPolicy {
                using Kernel = fem::mechanical::LinearElasticKernel<Element>;
                using DofTraits = fem::mechanical::LinearElasticDofTraits<Element::Dim>;
                using Material = material::LinearElasticMaterial<Element::Dim>;
                using MaterialTensor = typename Material::MaterialTensor;

                /// @brief Spatial dimension (2 or 3).
                static constexpr int Dim = Element::Dim;

                /**
                 * @brief Number of macroscopic strain load cases.
                 *
                 * The load cases are ordered in Voigt notation:
                 *
                 * - 2D: ε̄₁₁, ε̄₂₂, ε̄₁₂.
                 *
                 * - 3D: ε̄₁₁, ε̄₂₂, ε̄₃₃, ε̄₁₂, ε̄₁₃, ε̄₂₃.
                 */
                static constexpr int NumLoadCases = Material::VoigtSize;

                /// @brief Dof-major matrix type storing one global field per load case.
                using DofFieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, NumLoadCases>;
                
                /// @brief Nodal field type.
                using NodalFieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, Dim>;

                /// @brief List of nodal fields, one per load case.
                using NodalFieldList = std::array<NodalFieldMatrix, NumLoadCases>;

                /**
                 * @brief Linear elastic homogenization results.
                 *
                 * The homogenized stiffness tensor is always returned.
                 * Nodal displacements are stored only when requested
                 * through `FieldSave`.
                 */
                struct Results {
                    /// @brief Homogenized stiffness tensor.
                    MaterialTensor CBar;

                    /// @brief Total nodal displacements.
                    NodalFieldList u;

                    /// @brief Macroscopic nodal displacements.
                    NodalFieldList uMacro;

                    /// @brief Microscopic nodal displacements.
                    NodalFieldList uMicro;
                };

                /**
                 * @brief Macroscopic nodal displacements.
                 *
                 * For each prescribed macroscopic strain ε̄,
                 * the macroscopic displacements are:
                 *
                 * ```text
                 * ū=ε̄x
                 * ```
                 *
                 * @tparam Grid Grid type (e.g. Quad4Grid).
                 *
                 * @param[in] grid Grid.
                 *
                 * @returns Macroscopic nodal displacements.
                 */
                template <class Grid>
                static DofFieldMatrix makeMacroscopicFields(const Grid &grid) noexcept {
                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = DofTraits::NumNodeDofs * numNodes;

                    DofFieldMatrix UBar = DofFieldMatrix::Zero(static_cast<int>(numDofs), NumLoadCases);

                    for (std::size_t i = 0; i < numNodes; ++i) {
                        const auto node = grid.node(i);

                        if constexpr (Dim == 2) {
                            const double x = node(0);
                            const double y = node(1);

                            const int u = static_cast<int>(2 * i);
                            const int v = static_cast<int>(2 * i + 1);

                            // ū⁽¹¹⁾=ε̄⁽¹¹⁾x
                            UBar(u, 0) = x;
                            // ū⁽²²⁾=ε̄⁽²²⁾x
                            UBar(v, 1) = y;
                            // ū⁽¹²⁾=ε̄⁽¹²⁾x
                            UBar(u, 2) = 0.5 * y;
                            UBar(v, 2) = 0.5 * x;
                        }
                        else {
                            const double x = node(0);
                            const double y = node(1);
                            const double z = node(2);

                            const int u = static_cast<int>(3 * i);
                            const int v = static_cast<int>(3 * i + 1);
                            const int w = static_cast<int>(3 * i + 2);

                            // ū⁽¹¹⁾=ε̄⁽¹¹⁾x
                            UBar(u, 0) = x;
                            // ū⁽²²⁾=ε̄⁽²²⁾x
                            UBar(v, 1) = y;
                            // ū⁽³³⁾=ε̄⁽³³⁾x
                            UBar(w, 2) = z;
                            // ū⁽¹²⁾=ε̄⁽¹²⁾x
                            UBar(u, 3) = 0.5 * y;
                            UBar(v, 3) = 0.5 * x;
                            // ū⁽¹³⁾=ε̄⁽¹³⁾x
                            UBar(u, 4) = 0.5 * z;
                            UBar(w, 4) = 0.5 * x;
                            // ū⁽²³⁾=ε̄⁽²³⁾x
                            UBar(v, 5) = 0.5 * z;
                            UBar(w, 5) = 0.5 * y;
                        }
                    }

                    return UBar;
                }

                /**
                 * @brief Packages solver outputs into linear elastic results.
                 *
                 * @param[in] CBar Homogenized stiffness tensor.
                 * @param[in] U Total nodal displacements.
                 * @param[in] UMacro Macroscopic nodal displacements.
                 * @param[in] UMicro Microscopic nodal displacements.
                 * @param[in] fields Controls which nodal fields are stored.
                 *
                 * @returns Homogenization results.
                 */
                static Results makeResults(const MaterialTensor &CBar, const DofFieldMatrix &U, const DofFieldMatrix &UMacro, const DofFieldMatrix &UMicro, FieldSave fields) noexcept {
                    Results results;
                    results.CBar = CBar;

                    if (wants(fields, FieldSave::Total)) {
                        results.u = makeStoredFields_(U);
                    }

                    if (wants(fields, FieldSave::Macro)) {
                        results.uMacro = makeStoredFields_(UMacro);
                    }

                    if (wants(fields, FieldSave::Micro)) {
                        results.uMicro = makeStoredFields_(UMicro);
                    }

                    return results;
                }

            private:
                /**
                 * @brief Converts a nodal displacement block matrix into a list of nodal displacements.
                 *
                 * @param[in] U Nodal displacement block matrix.
                 *
                 * @return List of nodal displacements.
                 */
                static NodalFieldList makeStoredFields_(const DofFieldMatrix &U) {
                    const int numNodes = static_cast<int>(U.rows()) / Dim;
                    NodalFieldList out;

                    for (std::size_t i = 0; i < NumLoadCases; ++i) {
                        const NodalFieldMatrix u = U.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                        out[i] = u;
                    }

                    return out;
                }
            };

        } // namespace mechanical

    } // namespace solver

} // namespace monad
