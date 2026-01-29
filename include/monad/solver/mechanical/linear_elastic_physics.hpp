#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/fem/operator/mechanical/linear_elastic_matrix_free_operator_traits.hpp"
#include "monad/solver/solver_options.hpp"
#include "monad/solver/mechanical/linear_elastic_physics_traits.hpp"

namespace monad {

    namespace solver {

        namespace mechanical {

            /**
             * @brief Physics policy for linear elastic homogenization.
             *
             * @tparam Element Element class (e.g. Quad4).
             */
            template <class Element>
            struct LinearElasticPhysics {
                using Kernel = fem::mechanical::LinearElasticKernel<Element>;
                using Material = typename Kernel::Material;

                /// @brief Spatial dimension (2 or 3).
                static constexpr int Dim = Element::Dim;

                /**
                 * @brief Number of macroscopic fields:
                 *
                 * - D=2: ε̄₁₁, ε̄₂₂, ε̄₁₂.
                 *
                 * - D=3: ε̄₁₁, ε̄₂₂, ε̄₃₃, ε̄₁₂, ε̄₁₃, ε̄₂₃.
                 */
                static constexpr int NumMacroFields = Material::VoigtSize;

                using OperatorTraits = fem::mechanical::LinearElasticMatrixFreeOperatorTraits<Dim>;

                using Traits = LinearElasticPhysicsTraits<Dim>;

                /// @brief Global force matrix type.
                using FieldMatrix = Eigen::Matrix<double, Eigen::Dynamic, NumMacroFields>;

                /// @brief Homogenized stiffness tensor type.
                using MaterialTensor = typename Material::MaterialTensor;

                /**
                 * @brief Results from linear elastic homogenization.
                 *
                 * @note Displacement fields are stored in arrays
                 * ordered by the prescribed macroscopic strain
                 * loading directions:
                 *
                 * - D=2: ε̄₁₁, ε̄₂₂, ε̄₁₂.
                 *
                 * - D=3: ε̄₁₁, ε̄₂₂, ε̄₃₃, ε̄₁₂, ε̄₁₃, ε̄₂₃.
                 *
                 * Each array entry corresponds to one independent macroscopic loading case.
                 */
                struct Results {
                    using NodalField = Eigen::Matrix<double, Eigen::Dynamic, Dim>;

                    /// @brief Homogenized stiffness tensor C̄.
                    MaterialTensor CBar;

                    /// @brief Total nodal displacement fields u=ū+ũ.
                    std::array<NodalField, NumMacroFields> u;

                    /// @brief Macroscopic nodal displacement fields ū.
                    std::array<NodalField, NumMacroFields> uMacro;

                    /// @brief Microscopic nodal displacement fields ũ.
                    std::array<NodalField, NumMacroFields> uMicro;
                };

                /**
                 * @brief Macroscopic nodal displacements.
                 *
                 * The macroscopic displacements Ū are the displacements
                 * resulting from the macroscopic strains ε̄:
                 *
                 * Ū=ε̄x
                 *
                 * @tparam Grid Grid class (e.g. Quad4Grid).
                 *
                 * @param[in] grid Periodic unit cell grid.
                 *
                 * @returns Macroscopic nodal displacements.
                 */
                template <class Grid>
                static FieldMatrix macroscopicField(const GridBase<Grid, Element> &grid) {
                    return Traits::macroDisplacements(grid);
                }

                /**
                 * @brief Constructs linear elastic homogenization results.
                 *
                 * @param[in] CBar Homogenized stiffness tensor.
                 * @param[in] U Total nodal displacement fields.
                 * @param[in] UMacro Macroscopic nodal displacement fields.
                 * @param[in] UMicro Microscopic nodal displacement fields.
                 * @param[in] fields Controls which nodal fields are stored in the solver results.
                 *
                 * @returns Linear elastic homogenization results.
                 */
                static Results makeResults(const MaterialTensor &CBar, const FieldMatrix &U, const FieldMatrix &UMacro, const FieldMatrix &UMicro, FieldSave fields) noexcept {
                    using NodalField = typename Results::NodalField;

                    Results results;

                    results.CBar = CBar;

                    const auto numNodes = U.rows() / Dim;

                    if (wants(fields, FieldSave::Total)) {
                        std::array<NodalField, NumMacroFields> temp;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            temp[i] = U.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                        }

                        results.u = temp;
                    }

                    if (wants(fields, FieldSave::Macro)) {
                        std::array<NodalField, NumMacroFields> temp;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            temp[i] = UMacro.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                        }

                        results.uMacro = temp;
                    }

                    if (wants(fields, FieldSave::Micro)) {
                        std::array<NodalField, NumMacroFields> temp;

                        for (std::size_t i = 0; i < NumMacroFields; ++i) {
                            temp[i] = UMicro.col(static_cast<int>(i)).template reshaped<Eigen::RowMajor>(numNodes, Dim);
                        }

                        results.uMicro = temp;
                    }

                    return results;
                }
            };

        } // namespace mechanical

    } // namespace solver

} // namespace monad
