#pragma once

#include <array>
#include <cstddef>

namespace monad {

    namespace fem {

        namespace mechanical {

            /**
             * @brief Dof numbering rules for periodic linear elasticity problems.
             *
             * @tparam D Spatial dimension (2 or 3).
             */
            template <int D>
            struct LinearElasticDofTraits {
                /// @brief Number of dofs per node.
                static constexpr int NumNodeDofs = D;

                /// @brief Number of fixed periodic dofs to remove rigid-body transformations.
                static constexpr int NumFixedDofs = D;

                /**
                 * @brief Expands element node indices into dofs.
                 *
                 * Given an element with node indices [n₁, ..., nₖ], the
                 * associated dofs are
                 *
                 * - 2D: [2n₁, 2n₁+1, ..., 2nₖ, 2nₖ+1]
                 *
                 * - 3D: [3n₁, 3n₁+1, 3n₁+2, ..., 3nₖ, 3nₖ+1, 3nₖ+2]
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

                    std::array<std::size_t, NumNodeDofs * K> dofs;

                    std::size_t index = 0;

                    for (std::size_t i = 0; i < K; ++i) {
                        const std::size_t node = element[i];

                        for (std::size_t j = 0; j < NumNodeDofs; ++j) {
                            dofs[index++] = NumNodeDofs * node + j;
                        }
                    }

                    return dofs;
                }

                /**
                 * @brief Returns `true` if a periodic dof is fixed.
                 *
                 * The first `D` periodic dofs are fixed.
                 *
                 * - 2D: [0, 1]
                 *
                 * - 3D: [0, 1, 2]
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

                    return dof < NumFixedDofs;
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

        } // namespace mechanical

    } // namespace fem

} // namespace monad
