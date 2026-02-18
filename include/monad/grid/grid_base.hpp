#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <random>
#include <functional>
#include "monad/detail/constants.hpp"
#include "monad/integration/integrate_scalar.hpp"

namespace monad {

    /**
     * @brief Curiously recurring template pattern (CRTP) base class for grid meshes.
     *
     * @tparam Derived Concrete grid mesh class (CRTP).
     * @tparam ElementT Concrete element class (e.g. Quad4).
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
     *
     * - `void translate(const Resolution &shift) noexcept`
     */
    template <class Derived, class ElementT>
    class GridBase {
    public:
        static_assert(ElementT::Dim == 2 || ElementT::Dim == 3, "Element spatial dimension must be 2 or 3.");

        using Element = ElementT;

        /// @brief Spatial dimension (2 or 3).
        static constexpr int Dim = Element::Dim;

        using Resolution = std::array<std::size_t, Dim>;
        using Size = std::array<double, Dim>;

        using Point = typename Element::Point;
        using NodesMatrix = typename Element::NodesMatrix;

        using NodesList = std::vector<Point>;
        using ElementList = std::array<std::size_t, Element::NumNodes>;
        using ElementsList = std::vector<ElementList>;
        using DensityList = std::vector<double>;

        /**
         * @brief Constructs a grid mesh.
         *
         * @param[in] resolution How many cells in each dimension.
         * @param[in] size Physical lengths in each dimension.
         *
         * @throws std::invalid_argument if any entry in `resolution` is zero.
         * @throws std::invalid_argument if any entry in `size` is non-positive.
         *
         * @note Densities are initalized to 0.
         */
        GridBase(const Resolution &resolution, const Size &size)
            : resolution_(resolution), size_(size) {
            for (std::size_t i = 0; i < Dim; ++i) {
                if (resolution_[i] == 0) {
                    throw std::invalid_argument("Resolution in dimension " + std::to_string(i + 1) + " must be positive.");
                }

                if (size_[i] <= 0) {
                    throw std::invalid_argument("Size in dimension " + std::to_string(i + 1) + " (" + std::to_string(size_[i]) + ") must be positive.");
                }
            }

            densities_.resize(numElements());
            setDensitiesZeros();
        }

        /// @brief How many cells in each dimension.
        const Resolution &resolution() const noexcept {
            return resolution_;
        }

        /// @brief Physical lengths in each dimension.
        const Size &size() const noexcept {
            return size_;
        }

        /**
         * @brief Material density at each element.
         *
         * @note Densities are stored in row-major order.
         */
        const DensityList &densities() const noexcept {
            return densities_;
        }

        /// @brief Number of elements.
        std::size_t numElements() const noexcept {
            std::size_t product = 1;

            for (std::size_t entry : resolution_) {
                product *= entry;
            }

            return product;
        }

        /// @brief Number of nodes.
        std::size_t numNodes() const noexcept {
            return derived().numNodes();
        }

        /// @brief Number of periodic nodes.
        std::size_t numPeriodicNodes() const noexcept {
            return derived().numPeriodicNodes();
        }

        /**
         * @brief Coordinates for a specific node.
         *
         * @param[in] index Node index.
         *
         * @returns Coordinates.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numNodes()`).
         */
        Point node(std::size_t index) const {
            return derived().node(index);
        }

        /// @brief Coordinates for all nodes.
        NodesList nodes() const noexcept {
            NodesList data;
            data.reserve(numNodes());

            for (std::size_t i = 0; i < numNodes(); ++i) {
                data.push_back(node(i));
            }

            return data;
        }

        /**
         * @brief Node indices for a specific element.
         *
         * @param[in] index Element index.
         *
         * @returns Node indices.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         */
        ElementList element(std::size_t index) const {
            return derived().element(index);
        }

        /// @brief Node indices for all elements.
        ElementsList elements() const noexcept {
            ElementsList data;
            data.reserve(numElements());

            for (std::size_t i = 0; i < numElements(); ++i) {
                data.push_back(element(i));
            }

            return data;
        }

        /**
         * @brief Periodic node indices for a specific element.
         *
         * @param[in] index Element index.
         *
         * @returns Periodic node indices.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         */
        ElementList periodicElement(std::size_t index) const {
            return derived().periodicElement(index);
        }

        /// @brief Periodic node indices for all elements.
        ElementsList periodicElements() const noexcept {
            ElementsList data;
            data.reserve(numElements());

            for (std::size_t i = 0; i < numElements(); ++i) {
                data.push_back(periodicElement(i));
            }

            return data;
        }

        /**
         * @brief Material density for a specific element.
         *
         * @param[in] index Element index.
         *
         * @returns Material density.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         */
        double getDensity(std::size_t index) const {
            if (index >= numElements()) {
                throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
            }

            return densities_[index];
        }

        /**
         * @brief Set the material density for a specific element.
         *
         * @param[in] index Element index.
         * @param[in] density Element density.
         *
         * @throws std::invalid_argument if `density` is outside the range [0,1].
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         *
         * @note Densities are clamped to a minimum value `NUMERICAL_ZERO` to ensure numerical stability.
         */
        void setDensity(std::size_t index, double density) {
            if (index >= numElements()) {
                throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
            }

            if (density < 0 || density > 1) {
                throw std::invalid_argument("Density (" + std::to_string(density) + ") is out of range [0,1].");
            }

            densities_[index] = std::max(NUMERICAL_ZERO, density);
        }

        /**
         * @brief Set the material densities.
         *
         * @param[in] densities Material densities.
         *
         * @throws std::invalid_argument if the size of `densities` does not equal the number of grid elements.
         * @throws std::invalid_argument if any value in `densities` is outside the range [0,1].
         *
         * @note Densities are stored in row-major order.
         * @note Densities are clamped to a minimum value `NUMERICAL_ZERO` to ensure numerical stability.
         */
        void setDensities(DensityList densities) {
            if (densities.size() != numElements()) {
                throw std::invalid_argument("Densities size (" + std::to_string(densities.size()) + ") must equal number of grid elements (" + std::to_string(numElements()) + ").");
            }

            for (std::size_t i = 0; i < densities.size(); ++i) {
                setDensity(i, densities[i]);
            }
        }

        /**
         * @brief Set the material densities to a specific value.
         *
         * @param[in] density Element density.
         *
         * @throws std::invalid_argument if `density` is outside the range [0,1].
         *
         * @note Densities are clamped to a minimum value `NUMERICAL_ZERO` to ensure numerical stability.
         */
        void setDensitiesConstant(double density) {
            if (density < 0 || density > 1) {
                throw std::invalid_argument("Density (" + std::to_string(density) + ") is out of range [0,1].");
            }

            for (std::size_t i = 0; i < numElements(); ++i) {
                setDensity(i, density);
            }
        }

        /// @brief Set the material densities to zeros.
        void setDensitiesZeros() noexcept {
            setDensitiesConstant(0.0);
        }

        /// @brief Set the material densities to ones.
        void setDensitiesOnes() noexcept {
            setDensitiesConstant(1.0);
        }

        /**
         * @brief Set the material densities to random values.
         *
         * @param[in] seed RNG seed (default: -1 for random).
         *
         * @note Densities are clamped to a minimum value `NUMERICAL_ZERO` to ensure numerical stability.
         */
        void setDensitiesRandom(int seed = -1) noexcept {
            std::mt19937 rng;
            if (seed >= 0) {
                rng.seed(static_cast<unsigned int>(seed));
            }

            // std::uniform_real_distribution generates values in [0, 1), so setting
            // Use the next representable double after 1 to make the range effectively inclusive
            const double upperBound = std::nextafter(1.0, 2.0);

            std::uniform_real_distribution<double> dist(NUMERICAL_ZERO, upperBound);

            for (std::size_t i = 0; i < numElements(); ++i) {
                setDensity(i, dist(rng));
            }
        }

        /**
         * @brief Set the material densities using a continuous function.
         *
         * Element values are computed as the element-wise average of the continuous function.
         *
         * @param[in] f Continuous function from ℝᵈ→[0,1] evaluated at element nodes.
         *
         * @throws std::invalid_argument if `f` is outside the range [0,1].
         *
         * @note Densities are clamped to a minimum value `NUMERICAL_ZERO` to ensure numerical stability.
         */
        void setDensitiesFunction(const std::function<double(const Point &)> &f) {
            for (std::size_t i = 0; i < numElements(); ++i) {
                const auto nodes = elementNodes(i);

                auto integrand = [&](const Point &point) -> double {
                    const auto N = Element::shapeFunctions(point);
                    const Point globalPoint = N.transpose() * nodes;

                    const double density = f(globalPoint);

                    if (density < 0 || density > 1) {
                        throw std::invalid_argument("Function value (" + std::to_string(density) + ") is outside range [0,1].");
                    }

                    const auto J = Element::jacobian(point, nodes);

                    return density * std::abs(J.determinant());
                };

                double density = integration::integrateScalar(integrand, Element::quadratureRule()) / Element::measure(nodes);

                setDensity(i, density);
            }
        }

        /**
         * @brief Translates the element densities periodically.
         *
         * @param[in] shift Shift in each directions.
         */
        void translate(const Resolution &shift) noexcept {
            derived().translate(shift);
        }

        /**
         * @brief Nodal coordinates for a specific element.
         *
         * @param[in] index Element index.
         *
         * @returns Nodal coordinates.
         *
         * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
         */
        NodesMatrix elementNodes(std::size_t index) const {
            if (index >= numElements()) {
                throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
            }

            NodesMatrix data;

            int j = 0;
            for (std::size_t i : element(index)) {
                data.row(j++) = node(i);
            }

            return data;
        }

        /// @brief Grid area (2D) or volume (3D).
        double measure() const noexcept {
            const auto referenceElementNodes = elementNodes(0);
            const double referenceElementMeasure = Element::measure(referenceElementNodes);

            return referenceElementMeasure * static_cast<double>(numElements());
        }

        /// @brief Equality comparison.
        bool operator==(const GridBase<Derived, Element> &other) const noexcept {
            return resolution_ == other.resolution_
                && size_ == other.size_
                && densities_ == other.densities_;
        }

        /// @brief Inequality comparison.
        bool operator!=(const GridBase<Derived, Element> &other) const noexcept {
            return !(*this == other);
        }

    protected:
        /// @brief How many cells in each dimension.
        Resolution resolution_;

        /// @brief Physical lengths in each dimension.
        Size size_;

        /**
         * @brief Material density at each element.
         *
         * @note Densities are stored in row-major order.
         */
        DensityList densities_;

    private:
        /// @brief Utility method to safely cast the base class pointer to the derived class reference (CRTP).
        const Derived &derived() const noexcept {
            // The static_cast is safe because Derived must inherit from GridBase<Derived, Element>
            return *static_cast<const Derived *>(this);
        }
    };

} // namespace monad
