#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "monad/detail/constants.hpp"

namespace monad {

    namespace field {

        /**
         * @brief Per-element density field.
         *
         * @tparam D Spatial dimension (2 or 3).
         */
        template <int D>
        class DensityField {
        public:
            static_assert(D == 2 || D == 3, "Density field spatial dimension must be 2 or 3.");

            static constexpr int Dim = D;

            using Resolution = std::array<std::size_t, D>;
            using DensityList = std::vector<double>;

            /**
             * @brief Constructs a density field.
             *
             * @param[in] resolution Number of elements in each dimension.
             *
             * @throws std::invalid_argument if any entry in `resolution` is zero.
             *
             * @note Densities are initalized to zero.
             */
            explicit DensityField(const Resolution& resolution)
                : resolution_(resolution) {
                for (std::size_t i = 0; i < D; ++i) {
                    if (resolution_[i] == 0) {
                        throw std::invalid_argument("Resolution in dimension " + std::to_string(i + 1) + " must be positive.");
                    }
                }

                densities_.resize(numElements());
                setZeros();
            }

            /// @brief Number of elements in each dimension.
            const Resolution& resolution() const noexcept {
                return resolution_;
            }

            /// @brief Number of elements.
            std::size_t numElements() const noexcept {
                std::size_t total = 1;

                for (std::size_t n : resolution_) {
                    total *= n;
                }

                return total;
            }

            /**
             * @brief Per-element material densities.
             *
             * @note Densities are stored in row-major order.
             */
            const DensityList& densities() const noexcept {
                return densities_;
            }

            /**
             * @brief Material density of an element.
             *
             * @param[in] index Element index.
             *
             * @returns Material density of an element.
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
             * @brief Set the material density of an element.
             *
             * @param[in] index Element index.
             * @param[in] density Density.
             *
             * @throws std::invalid_argument if `density` is outside the range [0,1].
             * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
             *
             * @note Densities are clamped to a minimum value of `ZERO_TOLERANCE` to improve
             * numerical stability.
             */
            void setDensity(std::size_t index, double density) {
                if (index >= numElements()) {
                    throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
                }

                if (density < 0.0 || density > 1.0) {
                    throw std::invalid_argument("Density (" + std::to_string(density) + ") is out of range [0,1].");
                }

                densities_[index] = std::max(NUMERICAL_ZERO, density);
            }

            /**
             * @brief Set the material densities.
             *
             * @param[in] densities Densities.
             *
             * @throws std::invalid_argument if the size of `densities` does not equal the number of grid elements.
             * @throws std::invalid_argument if any value in `densities` is outside the range [0,1].
             *
             * @note Densities are stored in row-major order.
             * @note Densities are clamped to a minimum value of `ZERO_TOLERANCE` to improve numerical stability.
             */
            void setDensities(const DensityList& densities) {
                if (densities.size() != numElements()) {
                    throw std::invalid_argument("Densities size (" + std::to_string(densities.size()) + ") must equal number of elements (" + std::to_string(numElements()) + ").");
                }

                for (std::size_t i = 0; i < densities.size(); ++i) {
                    setDensity(i, densities[i]);
                }
            }

            /**
             * @brief Sets all material densities to a constant value.
             *
             * @param[in] density Density.
             *
             * @throws std::invalid_argument if `density` is outside the range [0,1].
             *
             * @note Densities are clamped to a minimum value of `ZERO_TOLERANCE` to
             * improve numerical stability.
             */
            void setConstant(double density) {
                if (density < 0.0 || density > 1.0) {
                    throw std::invalid_argument("Density (" + std::to_string(density) + ") is out of range [0,1].");
                }

                for (std::size_t i = 0; i < numElements(); ++i) {
                    setDensity(i, density);
                }
            }

            /// @brief Sets all material densities to zero.
            void setZeros() noexcept {
                setConstant(0.0);
            }

            /// @brief Sets all material densities to one.
            void setOnes() noexcept {
                setConstant(1.0);
            }

            /**
             * @brief Set the material densities to random values.
             *
             * @param[in] seed RNG seed (optional).
             *
             * @note Densities are clamped to a minimum value of
             * `ZERO_TOLERANCE` to improve numerical stability.
             */
            void setRandom(unsigned int seed = std::random_device{}()) noexcept {
                std::mt19937 rng(seed);

                for (std::size_t i = 0; i < numElements(); ++i) {
                    // std::uniform_real_distribution is not portably defined across compilers.
                    // Use the standardized mt19937 engine to guarantee consistency across Mac,
                    // Windows, and Linux.
                    const double density = static_cast<double>(rng()) / static_cast<double>(std::mt19937::max());

                    setDensity(i, density);
                }
            }

            /**
             * @brief Periodically translates the element densities.
             *
             * @param[in] shift Shift in each directions.
             */
            void translate(const Resolution& shift) noexcept {
                DensityList shifted(numElements());

                if constexpr (D == 2) {
                    const std::size_t nx = resolution_[0];
                    const std::size_t ny = resolution_[1];

                    for (std::size_t i = 0; i < nx; ++i) {
                        for (std::size_t j = 0; j < ny; ++j) {
                            const std::size_t oldIndex = j * nx + i;
                            const std::size_t iNew = (i + shift[0]) % nx;
                            const std::size_t jNew = (j + shift[1]) % ny;
                            const std::size_t newIndex = jNew * nx + iNew;
                            
                            shifted[newIndex] = densities_[oldIndex];
                        }
                    }
                }
                else {
                    const std::size_t nx = resolution_[0];
                    const std::size_t ny = resolution_[1];
                    const std::size_t nz = resolution_[2];
                    const std::size_t elementsPerPlane = nx * ny;

                    for (std::size_t i = 0; i < nx; ++i) {
                        for (std::size_t j = 0; j < ny; ++j) {
                            for (std::size_t k = 0; k < nz; ++k) {
                                const std::size_t oldIndex = k * elementsPerPlane + j * nx + i;
                                const std::size_t iNew = (i + shift[0]) % nx;
                                const std::size_t jNew = (j + shift[1]) % ny;
                                const std::size_t kNew = (k + shift[2]) % nz;
                                const std::size_t newIndex = kNew * elementsPerPlane + jNew * nx + iNew;
                                
                                shifted[newIndex] = densities_[oldIndex];
                            }
                        }
                    }
                }

                setDensities(shifted);
            }

            /// @brief Equality comparison.
            bool operator==(const DensityField& other) const noexcept {
                auto vectorEqual = [](const DensityList& v1, const DensityList& v2) noexcept {
                    if (v1.size() != v2.size()) {
                        return false;
                    }

                    for (std::size_t i = 0; i < v1.size(); ++i) {
                        if (std::fabs(v1[i] - v2[i]) > NUMERICAL_ZERO) {
                            return false;
                        }
                    }

                    return true;
                };

                return resolution_ == other.resolution_ && vectorEqual(densities_, other.densities_);
            }

            /// @brief Inequality comparison.
            bool operator!=(const DensityField& other) const noexcept {
                return !(*this == other);
            }
            
        private:
            /// @brief Number of elements in each dimension.
            Resolution resolution_;

            /// @brief Per-element material densities.
            DensityList densities_;
        };

    } // namespace field

} // namespace monad
