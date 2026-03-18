#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

namespace monad {

    namespace material {

        /**
         * @brief Linear elastic material model.
         *
         * By Hooke's law, the stress σ∈ℝᵛ is a linear function of
         * the strain ε∈ℝᵛ:
         *
         * ```text
         * σ=Cε
         * ```
         *
         * - C∈Symᵥ(ℝ) is the stiffness tensor in Voigt notation.
         *
         * @tparam D Spatial dimension (2 or 3).
         */
        template <int D>
        class LinearElasticMaterial {
        public:
            static_assert(D == 2 || D == 3, "Spatial dimension D must be 2 or 3.");

            /// @brief Spatial dimension (2 or 3).
            static constexpr int Dim = D;

            /// @brief Number of components in Voigt notation.
            static constexpr int VoigtSize = (Dim == 2) ? 3 : 6;

            /// @brief Stiffness tensor type.
            using MaterialTensor = Eigen::Matrix<double, VoigtSize, VoigtSize>;

            /// @brief Default constructor.
            LinearElasticMaterial() = default;

            /**
             * @brief Constructs a linear elastic material.
             *
             * @param[in] C Stiffness tensor.
             *
             * @throws std::invalid_argument if `C` is not positive definite.
             */
            explicit LinearElasticMaterial(const MaterialTensor &C)
                : C_(C) {
                if (!detail::isPD(C_)) {
                    throw std::invalid_argument("Stiffness tensor is not positive definite.");
                }
            }

            /// @brief Stiffness tensor.
            const MaterialTensor &materialTensor() const noexcept {
                return C_;
            }

            /// @brief Equality comparison.
            bool operator==(const LinearElasticMaterial &other) const noexcept {
                return C_.isApprox(other.C_, NUMERICAL_ZERO);
            }

            /// @brief Inequality comparison.
            bool operator!=(const LinearElasticMaterial &other) const noexcept {
                return !(*this == other);
            }

        protected:
            /// @brief Stiffness tensor.
            MaterialTensor C_;
        };

    } // namespace material

} // namespace monad
