#pragma once

#include <array>
#include <cstddef>

namespace monad {

    namespace fem {

        namespace scalar {

            /// @brief Dof numbering rules for periodic linear scalar diffusive problems.
            struct LinearScalarDiffusiveDofTraits {
                /// @brief Number of dofs per node.
                static constexpr int NumNodeDofs = 1;

                /// @brief Number of fixed periodic dofs to remove rigid-body transformations.
                static constexpr int NumFixedDofs = 1;

                /**
                 * @brief Expands element node indices into dofs.
                 *
                 * Given an element with node indices [n₁, ..., nₖ], the
                 * associated dofs are [n₁, ..., nₖ].
                 *
                 * @tparam K Number of element nodes.
                 *
                 * @param[in] element Element node indices.
                 * @param[in] numNodes Total number of nodes.
                 *
                 * @returns Dofs associated with `element`.
                 *
                 * @note `numNodes` is unused here. It is kept only to match
                 * the interface of other dof traits that require node-based offsets.
                 */
                template <std::size_t K>
                static auto elementDofs(const std::array<std::size_t, K> &element, std::size_t numNodes) noexcept {
                    (void)numNodes;

                    return element;
                }

                /**
                 * @brief Returns `true` if a periodic dof is fixed.
                 *
                 * The first periodic dof is fixed.
                 *
                 * @param[in] dof Periodic dof.
                 * @param[in] numNodes Total number of periodic nodes.
                 *
                 * @returns `true` if `dof` is fixed, `false` otherwise.
                 *
                 * @note `numNodes` is unused here. It is kept only to match
                 * the interface of other dof traits that require node-based offsets.
                 */
                static bool isFixedPeriodicDof(std::size_t dof, std::size_t numNodes) noexcept {
                    (void)numNodes;

                    return dof == 0;
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
                 *
                 * @note `numNodes` is unused here. It is kept only to match
                 * the interface of other dof traits that require node-based offsets.
                 */
                static std::size_t periodicToReducedDof(std::size_t dof, std::size_t numNodes) noexcept {
                    (void)numNodes;

                    return dof - NumFixedDofs;
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
                 *
                 * @note `numNodes` is unused here. It is kept only to match
                 * the interface of other dof traits that require node-based offsets.
                 */
                static std::size_t reducedToPeriodicDof(std::size_t dof, std::size_t numNodes) noexcept {
                    (void)numNodes;

                    return dof + NumFixedDofs;
                }
            };

        } // namespace scalar

    } // namespace fem

} // namespace monad
