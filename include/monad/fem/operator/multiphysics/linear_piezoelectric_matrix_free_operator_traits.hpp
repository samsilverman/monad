#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/operator/mechanical/linear_elastic_matrix_free_operator_traits.hpp"
#include "monad/fem/operator/scalar/linear_scalar_diffusive_matrix_free_operator_traits.hpp"

namespace monad {

    namespace fem {

        namespace multiphysics {

            /**
             * @brief Linear piezoelectric matrix free operator traits.
             *
             * This struct isolates physics- and dimension-specific logic
             * so that the matrix free operator can be written generically.
             *
             * @tparam D Spatial dimension (2 or 3).
             */
            template <int D>
            struct LinearPiezoelectricMatrixFreeOperatorTraits {
                using MechanicalTraits = mechanical::LinearElasticMatrixFreeOperatorTraits<D>;
                using ElectricalTraits = scalar::LinearScalarDiffusiveMatrixFreeOperatorTraits;

                /// @brief Number of dofs per node.
                static constexpr int NumNodeDofs = MechanicalTraits::NumNodeDofs + ElectricalTraits::NumNodeDofs;

                /// @brief Number of fixed dofs to remove rigid body transformations.
                static constexpr int NumFixedDofs = MechanicalTraits::NumFixedDofs + ElectricalTraits::NumFixedDofs;

                /**
                 * @brief List of dofs associated with an element's nodes.
                 *
                 * Given an element with node indices [n₁, ..., nₖ], its dofs are:
                 *
                 * - D=2: [2n₁, 2n₁+1, ..., 2nₖ, 2nₖ+1, n₁+2k, ..., nₖ+2k]
                 *
                 * - D=3: [3n₁, 3n₁+1, 3n₁+2, ..., 3nₖ, 3nₖ+1, 3nₖ+2, n₁+3k, ..., nₖ+3k]
                 *
                 * @tparam K Number of nodes in the element.
                 *
                 * @param[in] element Node indices for a specific element.
                 * @param[in] numNodes Number of nodes.
                 *
                 * @returns List of dofs associated with `element`'s nodes.
                 *
                 * @note The `numNodes` argument is included for multi-physics
                 * operators where dof offsets may depend on the total number
                 * of nodes.
                 */
                template <std::size_t K>
                static auto dofs(const std::array<std::size_t, K> &element, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;

                    // Mechanical dofs
                    const auto mechanicalDofs = MechanicalTraits::dofs(element, numNodes);

                    // Electrical dofs
                    auto electricalDofs = ElectricalTraits::dofs(element, numNodes);

                    // Offset electrical dofs by number of mechanical dofs
                    for (auto &dof : electricalDofs) {
                        dof += numMechanicalDofs;
                    }

                    // Concatenate
                    std::array<std::size_t, NumNodeDofs * K> data;
                    std::size_t index = 0;

                    for (auto dof : mechanicalDofs) {
                        data[index++] = dof;
                    }
                    for (auto dof : electricalDofs) {
                        data[index++] = dof;
                    }

                    return data;
                }

                /**
                 * @brief Checks if a global dof is fixed.
                 *
                 * The fixed dofs are (assuming m total nodes):
                 *
                 * - D=2: [0, 1, 2m]
                 *
                 * - D=3: [0, 1, 2, 3m]
                 *
                 * @param[in] dof Global dof.
                 * @param[in] numNodes Number of nodes.
                 *
                 * @returns `true` if `dof` is fixed, `false` otherwise.
                 *
                 * @note The `numNodes` argument is included for multi-physics
                 * operators where dof offsets may depend on the total number
                 * of nodes.
                 */
                static bool isFixedDof(std::size_t dof, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;

                    // Fixed mechanical dofs
                    if (dof < numMechanicalDofs) {
                        return MechanicalTraits::isFixedDof(dof, numNodes);
                    }

                    // Fixed electrical dof
                    // Electrical dofs are offset by number of mechanical dofs
                    return ElectricalTraits::isFixedDof(dof - numMechanicalDofs, numNodes);
                }

                /**
                 * @brief Maps a global dof to its reduced version.
                 *
                 * Matrix-free operators act only on unconstrained periodic dofs.
                 * This function maps dofs from a space where fixed dofs are
                 * included to a space where they are removed.
                 *
                 * @param[in] dof Global dof.
                 * @param[in] numNodes Number of nodes.
                 *
                 * @returns Reduced dof.
                 *
                 * @note The `numNodes` argument is included for multi-physics
                 * operators where dof offsets may depend on the total number
                 * of nodes.
                 */
                static std::size_t reducedDof(std::size_t dof, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;
                    const std::size_t numReducedMechanicalDofs = numMechanicalDofs - MechanicalTraits::NumFixedDofs;

                    // Shifting for mechanical dofs.
                    if (dof < numMechanicalDofs) {
                        return MechanicalTraits::reducedDof(dof, numNodes);
                    }

                    // Shifting for electrical dof
                    // Electrical dofs are offset by number of mechanical dofs
                    return ElectricalTraits::reducedDof(dof - numMechanicalDofs, numNodes) + numReducedMechanicalDofs;
                }

                /**
                 * @brief Maps a reduced dof to its global version.
                 *
                 * Matrix-free operators act only on unconstrained periodic dofs.
                 * This function maps dofs from a space where fixed dofs
                 * are removed to a space where they are included.
                 *
                 * @param[in] dof Reduced dof.
                 * @param[in] numNodes Number of nodes.
                 *
                 * @returns Global dof.
                 *
                 * @note The `numNodes` argument is included for multi-physics
                 * operators where dof offsets may depend on the total number
                 * of nodes.
                 */
                static std::size_t expandedDof(std::size_t dof, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;
                    const std::size_t numReducedMechanicalDofs = numMechanicalDofs - MechanicalTraits::NumFixedDofs;

                    // Shifting for mechanical dofs
                    if (dof < numReducedMechanicalDofs) {
                        return MechanicalTraits::expandedDof(dof, numNodes);
                    }

                    // Shifting for electrical dof
                    // Electrical dofs are offset by number of mechanical dofs
                    return ElectricalTraits::expandedDof(dof - numReducedMechanicalDofs, numNodes) + numMechanicalDofs;
                }
            };

        } // namespace multiphysics

    } // namespace fem

} // namespace monad
