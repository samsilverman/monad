#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace monad {

    namespace fem {

        /**
         * @brief Maps element periodic dofs into the reduced periodic dof space.
         *
         * Fixed dofs are stored as `-1`.
         *
         * @tparam Grid Grid type (e.g. Quad4Grid).
         * @tparam Traits Dof traits type (e.g. LinearElasticDofTraits).
         */
        template <class Grid, class Traits>
        class DofMap {
        public:
            using Element = typename Grid::Element;
            using ElementDofs = std::array<int, Element::NumNodes * Traits::NumNodeDofs>;
            using ElementDofList = std::vector<ElementDofs>;

            /**
             * @brief Constructs a dof map.
             *
             * @param[in] grid Grid.
             */
            explicit DofMap(const Grid &grid) {
                numPeriodicNodes_ = grid.numPeriodicNodes();
                numReducedDofs_ = Traits::NumNodeDofs * numPeriodicNodes_ - Traits::NumFixedDofs;

                elementDofs_.reserve(grid.numElements());

                for (std::size_t i = 0; i < grid.numElements(); ++i) {
                    const auto element = grid.periodicElement(i);
                    const auto periodicDofs = Traits::elementDofs(element, numPeriodicNodes_);

                    ElementDofs reducedDofs;

                    for (std::size_t j = 0; j < periodicDofs.size(); ++j) {
                        const std::size_t periodicDof = periodicDofs[j];

                        if (Traits::isFixedPeriodicDof(periodicDof, numPeriodicNodes_)) {
                            reducedDofs[j] = -1;
                        }
                        else {
                            const std::size_t reducedDof = Traits::periodicToReducedDof(periodicDof, numPeriodicNodes_); 
                            reducedDofs[j] = static_cast<int>(reducedDof);
                        }
                    }

                    elementDofs_.push_back(reducedDofs);
                }
            }

            /// @brief Number of elements in the map.
            std::size_t size() const noexcept {
                return elementDofs_.size();
            }

            /// @brief Number of reduced periodic dofs.
            std::size_t numReducedDofs() const noexcept {
                return numReducedDofs_;
            }

            /**
             * @brief Reduced periodic dofs of an element.
             *
             * @param[in] index Element index.
             *
             * @returns Reduced periodic dofs of an element.
             *
             * @throws std::out_of_range if `index` is outside the range [0, `size()`).
             */
            const ElementDofs &reducedDofs(std::size_t index) const {
                if (index >= elementDofs_.size()) {
                    throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(elementDofs_.size()) + ").");
                }

                return elementDofs_[index];
            }

            /**
             * @brief Returns `true` if a periodic dof is fixed.
             *
             * @param[in] dof Periodic dof.
             *
             * @returns `true` if `dof` is fixed, `false` otherwise.
             */
            bool isFixedPeriodicDof(std::size_t dof) const noexcept {
                return Traits::isFixedPeriodicDof(dof, numPeriodicNodes_);
            }

            /**
             * @brief Maps a periodic dof to a reduced periodic dof.
             *
             * The reduced periodic dof space is formed by removing fixed
             * periodic dofs from the full periodic dof numbering.
             *
             * @param[in] dof Periodic dof.
             *
             * @returns Reduced periodic dof.
             */
            std::size_t periodicToReducedDof(std::size_t dof) const noexcept {
                return Traits::periodicToReducedDof(dof, numPeriodicNodes_);
            }

            /**
             * @brief Maps a reduced periodic dof to a periodic dof.
             *
             * The reduced periodic dof space is formed by removing fixed
             * periodic dofs from the full periodic dof numbering.
             *
             * @param[in] dof Reduced periodic dof.
             *
             * @returns Periodic dof.
             */
            std::size_t reducedToPeriodicDof(std::size_t dof) const noexcept {
                return Traits::reducedToPeriodicDof(dof, numPeriodicNodes_);
            }

        private:
            /// @brief Number of periodic nodes.
            std::size_t numPeriodicNodes_;

            /// @brief Number of reduced periodic dofs.
            std::size_t numReducedDofs_;

            /// @brief Per-element mapping into the reduced periodic dof space.
            ElementDofList elementDofs_;
        };

    } // namespace fem

} // namespace monad
