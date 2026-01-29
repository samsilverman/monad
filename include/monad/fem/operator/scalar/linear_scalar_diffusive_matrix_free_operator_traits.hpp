#pragma once

#include <array>
#include <cstddef>

namespace monad {

    namespace fem {

        namespace scalar {

            /**
             * @brief Linear scalar diffusive matrix free operator traits.
             *
             * This struct isolates physics- and dimension-specific logic
             * so that the matrix free operator can be written generically.
             */
            struct LinearScalarDiffusiveMatrixFreeOperatorTraits {
                /// @brief Number of dofs per node.
                static constexpr int NumNodeDofs = 1;

                /// @brief Number of fixed dofs to remove rigid body transformations.
                static constexpr int NumFixedDofs = 1;

                /**
                 * @brief List of dofs associated with an element's nodes.
                 *
                 * Given an element with node indices [n₁, ..., nₖ], its dofs are:
                 *
                 * [n₁, ..., nₖ]
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
                    (void)numNodes;

                    auto dofs = element;

                    return dofs;
                }

                /**
                 * @brief Checks if a global dof is fixed.
                 *
                 * The fixed dofs are:
                 *
                 * [0]
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
                    (void)numNodes;

                    return dof == 0;
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
                    (void)numNodes;

                    return dof - NumFixedDofs;
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
                    (void)numNodes;

                    return dof + NumFixedDofs;
                }
            };

        } // namespace scalar

    } // namespace fem

} // namespace monad
