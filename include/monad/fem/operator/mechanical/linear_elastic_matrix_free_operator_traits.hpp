#pragma once

#include <array>
#include <cstddef>

namespace monad {

    namespace fem {

        namespace mechanical {

            /**
             * @brief Linear elastic matrix free operator traits.
             *
             * This struct isolates physics- and dimension-specific logic
             * so that the matrix free operator can be written generically.
             *
             * @tparam D Spatial dimension (2 or 3).
             */
            template <int D>
            struct LinearElasticMatrixFreeOperatorTraits {
                /// @brief Number of dofs per node.
                static constexpr int NumNodeDofs = D;

                /// @brief Number of fixed dofs to remove rigid body transformations.
                static constexpr int NumFixedDofs = D;

                /**
                 * @brief List of dofs associated with an element's nodes.
                 *
                 * Given an element with node indices [n₁, ..., nₖ], its dofs are:
                 *
                 * - D=2: [2n₁, 2n₁+1, ..., 2nₖ, 2nₖ+1]
                 *
                 * - D=3: [3n₁, 3n₁+1, 3n₁+2, ..., 3nₖ, 3nₖ+1, 3nₖ+2]
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

                    std::array<std::size_t, NumNodeDofs * K> data;

                    std::size_t index = 0;

                    for (std::size_t i = 0; i < K; ++i) {
                        const std::size_t node = element[i];

                        for (std::size_t j = 0; j < NumNodeDofs; ++j) {
                            data[index++] = NumNodeDofs * node + j;
                        }
                    }

                    return data;
                }

                /**
                 * @brief Checks if a global dof is fixed.
                 *
                 * The fixed dofs are:
                 *
                 * - D=2: [0, 1]
                 *
                 * - D=3: [0, 1, 2]
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
                static bool isFixedDof(std::size_t dof, std::size_t numNodes) noexcept
                {
                    (void)numNodes;

                    return dof < NumFixedDofs;
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

        } // namespace mechanical

    } // namespace fem

} // namespace monad
