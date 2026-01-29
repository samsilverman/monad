#pragma once

#include "monad/grid/grid_base.hpp"

namespace monad {

    /**
     * @brief Curiously recurring template pattern (CRTP) base class for 3D grid meshes.
     *
     * @tparam Derived Concrete grid mesh class (CRTP).
     * @tparam Element Concrete element class (e.g. Hex8).
     *
     * @note `Derived` classes must provide:
     *
     * - `std::size_t numNodes() const noexcept`
     *
     * - `std::size_t numPeriodicNodes() const noexcept`
     *
     * - `Point node(std::size_t index) const`
     *
     * - `ElementList element(std::size_t index) const`
     *
     * - `ElementList periodicElement(std::size_t index) const`
     */
    template <class Derived, class Element>
    class Grid3dBase : public GridBase<Derived, Element> {
    public:
        static_assert(Element::Dim == 3, "Element spatial dimension must be 3.");

        using GridBase<Derived, Element>::GridBase;
        using Resolution = typename GridBase<Derived, Element>::Resolution;
        using DensityList = typename GridBase<Derived, Element>::DensityList;

        /**
         * @brief Translates the element densities periodically.
         *
         * @param[in] shift Shift in each directions.
         */
        void translate(const Resolution &shift) noexcept {
            DensityList shiftedDensities(this->numElements());

            const std::size_t nx = this->resolution_[0];
            const std::size_t ny = this->resolution_[1];
            const std::size_t nz = this->resolution_[2];

            const std::size_t elementsPerPlane = nx * ny;

            for (std::size_t i = 0; i < nx; ++i) {
                for (std::size_t j = 0; j < ny; ++j) {
                    for (std::size_t k = 0; k < nz; ++k) {
                        const std::size_t oldIndex = k * elementsPerPlane + j * nx + i;

                        const std::size_t iNew = (i + shift[0]) % nx;
                        const std::size_t jNew = (j + shift[1]) % ny;
                        const std::size_t kNew = (k + shift[2]) % nz;

                        std::size_t newIndex = kNew * elementsPerPlane + jNew * nx + iNew;

                        shiftedDensities[newIndex] = this->densities_[oldIndex];
                    }
                }
            }

            this->setDensities(shiftedDensities);
        }

        /// @brief Grid volume.
        double volume() const noexcept {
            return this->measure();
        }
    };

} // namespace monad
