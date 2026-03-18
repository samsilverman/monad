#pragma once

#include <stdexcept>
#include <string>
#include <Eigen/Core>
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

namespace monad {

    namespace material {

        /**
         * @brief Linear transport material model.
         *
         * In a linear transport constitutive law, the flux J∈ℝᵈ is
         * a linear function of a scalar potential gradient ∇φ∈ℝᵈ:
         *
         * ```text
         * J=-K∇φ
         * ```
         *
         * - K∈Symᴅ(ℝ) is the transport tensor.
         *
         * @tparam D Spatial dimension (2 or 3).
         */
        template <int D>
        class LinearTransportMaterial {
        public:
            static_assert(D == 2 || D == 3, "Spatial dimension D must be 2 or 3.");

            /// @brief Spatial dimension (2 or 3).
            static constexpr int Dim = D;

            /// @brief Transport tensor type.
            using MaterialTensor = Eigen::Matrix<double, Dim, Dim>;

            /**
             * @brief Constructs a linear transport material.
             *
             * @param[in] K Transport tensor.
             *
             * @throws std::invalid_argument if `K` is not positive definite.
             */
            explicit LinearTransportMaterial(const MaterialTensor &K)
                : K_(K) {
                if (!detail::isPD(K_)) {
                    throw std::invalid_argument("Transport tensor is not positive definite.");
                }
            }

            /**
             * @brief Constructs an isotropic linear transport material.
             *
             * @param[in] K Transport constant.
             *
             * @throws std::invalid_argument if `K` is non-positive.
             */
            explicit LinearTransportMaterial(double K) {
                if (K <= 0.0) {
                    throw std::invalid_argument("K (" + std::to_string(K) + ") must be positive.");
                }

                K_ = MaterialTensor::Identity() * K;
            }

            /// @brief Transport tensor.
            const MaterialTensor &materialTensor() const noexcept {
                return K_;
            }

            /// @brief Equality comparison.
            bool operator==(const LinearTransportMaterial &other) const noexcept {
                return K_.isApprox(other.K_, NUMERICAL_ZERO);
            }

            /// @brief Inequality comparison.
            bool operator!=(const LinearTransportMaterial &other) const noexcept {
                return !(*this == other);
            }

        private:
            /// @brief Transport tensor.
            MaterialTensor K_;
        };

    } // namespace material

} // namespace monad
