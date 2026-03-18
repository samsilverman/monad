#pragma once

#include <array>
#include <cstddef>
#include "monad/fem/operator/mechanical/linear_elastic_dof_traits.hpp"
#include "monad/fem/operator/scalar/linear_scalar_diffusive_dof_traits.hpp"

namespace monad {

    namespace fem {

        namespace multiphysics {

            /**
             * @brief Dof numbering rules for periodic linear piezoelectricity problems.
             *
             * @tparam D Spatial dimension (2 or 3).
             */
            template <int D>
            struct LinearPiezoelectricDofTraits {
                using MechanicalTraits = mechanical::LinearElasticDofTraits<D>;
                using ElectricalTraits = scalar::LinearScalarDiffusiveDofTraits;

                /// @brief Number of dofs per node.
                static constexpr int NumNodeDofs = MechanicalTraits::NumNodeDofs + ElectricalTraits::NumNodeDofs;

                /// @brief Number of fixed periodic dofs to remove rigid-body transformations.
                static constexpr int NumFixedDofs = MechanicalTraits::NumFixedDofs + ElectricalTraits::NumFixedDofs;

                /**
                 * @brief Expands element node indices into dofs.
                 *
                 * Given an element with node indices [n₁, ..., nₖ], the
                 * associated dofs are
                 *
                 * - 2D: [2n₁, 2n₁+1, ..., 2nₖ, 2nₖ+1, n₁+2k, ..., nₖ+2k]
                 *
                 * - 3D: [3n₁, 3n₁+1, 3n₁+2, ..., 3nₖ, 3nₖ+1, 3nₖ+2, n₁+3k, ..., nₖ+3k]
                 *
                 * @tparam K Number of element nodes.
                 *
                 * @param[in] element Element node indices.
                 * @param[in] numNodes Total number of nodes.
                 *
                 * @returns Dofs associated with `element`.
                 */
                template <std::size_t K>
                static auto elementDofs(const std::array<std::size_t, K> &element, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;

                    // Mechanical dofs
                    const auto mechanicalDofs = MechanicalTraits::elementDofs(element, numNodes);

                    // Electrical dofs
                    auto electricalDofs = ElectricalTraits::elementDofs(element, numNodes);

                    // Offset electrical dofs by number of mechanical dofs
                    for (auto &dof : electricalDofs) {
                        dof += numMechanicalDofs;
                    }

                    // Concatenate
                    std::array<std::size_t, NumNodeDofs * K> dofs;
                    std::size_t index = 0;

                    for (auto dof : mechanicalDofs) {
                        dofs[index++] = dof;
                    }
                    for (auto dof : electricalDofs) {
                        dofs[index++] = dof;
                    }

                    return dofs;
                }

                /**
                 * @brief Returns `true` if a periodic dof is fixed.
                 *
                 * The first `D` mechanical dofs and the first electrical
                 * dof are fixed:
                 *
                 * - 2D: [0, 1, 2m]
                 *
                 * - 3D: [0, 1, 2, 3m]
                 *
                 * @param[in] dof Periodic dof.
                 * @param[in] numNodes Total number of periodic nodes.
                 *
                 * @returns `true` if `dof` is fixed, `false` otherwise.
                 */
                static bool isFixedPeriodicDof(std::size_t dof, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;

                    // Fixed mechanical periodic dofs
                    if (dof < numMechanicalDofs) {
                        return MechanicalTraits::isFixedPeriodicDof(dof, numNodes);
                    }

                    // Fixed electrical periodic dof
                    // Electrical periodic dofs are offset by number of mechanical periodic dofs
                    return ElectricalTraits::isFixedPeriodicDof(dof - numMechanicalDofs, numNodes);
                }

                /**
                 * @brief Maps a periodic dof to a reduced periodic dof.
                 *
                 * The reduced periodic dof space is formed by removing fixed
                 * periodic dofs from the full periodic dof numbering.
                 *
                 * @param[in] dof Periodic dof.
                 * @param[in] numNodes Total number of periodic nodes.
                 *
                 * @returns Reduced periodic dof.
                 */
                static std::size_t periodicToReducedDof(std::size_t dof, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;
                    const std::size_t numReducedMechanicalDofs = numMechanicalDofs - MechanicalTraits::NumFixedDofs;

                    // Mechanical periodic dofs
                    if (dof < numMechanicalDofs) {
                        return MechanicalTraits::periodicToReducedDof(dof, numNodes);
                    }

                    // Electrical periodic dofs
                    // Electrical periodic dofs are offset by number of mechanical periodic dofs
                    return ElectricalTraits::periodicToReducedDof(dof - numMechanicalDofs, numNodes) + numReducedMechanicalDofs;
                }

                /**
                 * @brief Maps a reduced periodic dof to a periodic dof.
                 *
                 * The reduced periodic dof space is formed by removing fixed
                 * periodic dofs from the full periodic dof numbering.
                 *
                 * @param[in] dof Reduced periodic dof.
                 * @param[in] numNodes Total number of periodic nodes.
                 *
                 * @returns Periodic dof.
                 */
                static std::size_t reducedToPeriodicDof(std::size_t dof, std::size_t numNodes) noexcept {
                    const std::size_t numMechanicalDofs = numNodes * MechanicalTraits::NumNodeDofs;
                    const std::size_t numReducedMechanicalDofs = numMechanicalDofs - MechanicalTraits::NumFixedDofs;

                    // Mechanical reduced periodic dofs
                    if (dof < numReducedMechanicalDofs) {
                        return MechanicalTraits::reducedToPeriodicDof(dof, numNodes);
                    }

                    // Electrical reduced periodic dofs
                    // Electrical reduced periodic dofs are offset by number of mechanical reduced periodic dofs
                    return ElectricalTraits::reducedToPeriodicDof(dof - numReducedMechanicalDofs, numNodes) + numMechanicalDofs;
                }
            };

        } // namespace multiphysics

    } // namespace fem

} // namespace monad
