#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/detail/eigen_utils.hpp"

namespace monad {

    namespace material {

        /**
         * @brief Linear piezoelectric material model.
         *
         * In the stress-charge form, mechanical fields (stress T∈ℝᵛ and strain S∈ℝᵛ)
         * and electrical fields (electric displacement D∈ℝᵈ and electric field E∈ℝᵈ)
         * are coupled by:
         *
         * ```text
         * T=cS-dᵀE
         * D=dT+ϵE
         * ```
         *
         * - c∈Symᵥ(ℝ) is the stiffness tensor in Voigt notation.
         *
         * - ϵ∈Symᴅ(ℝ) is the permittivity tensor.
         *
         * - d∈ℝᵈˣᵛ is the piezoelectric coupling tensor.
         *
         * @tparam D Spatial dimension (2 or 3).
         */
        template <int D>
        class LinearPiezoelectricMaterial {
        public:
            static_assert(D == 2 || D == 3, "Spatial dimension D must be 2 or 3.");

            /// @brief Spatial dimension (2 or 3).
            static constexpr int Dim = D;

            /// @brief Number of components in Voigt notation.
            static constexpr int VoigtSize = (Dim == 2) ? 3 : 6;

            /// @brief Stiffness tensor type.
            using StiffnessTensor = Eigen::Matrix<double, VoigtSize, VoigtSize>;

            /// @brief Permittivity tensor type.
            using PermittivityTensor = Eigen::Matrix<double, Dim, Dim>;

            /// @brief Piezoelectric coupling tensor type.
            using CouplingTensor = Eigen::Matrix<double, Dim, VoigtSize>;

            /** @brief Coupled constitutive operator type.
             *
             * ```text
             * ⎡ c  -dᵀ⎤
             * ⎣-d  -ϵ ⎦
             * ```
             */
            using MaterialTensor = Eigen::Matrix<double, VoigtSize + Dim, VoigtSize + Dim>;

            /**
             * @brief Constructs a linear piezoelectric material.
             *
             * @param[in] c Stiffness tensor.
             * @param[in] epsilon Permittivity tensor.
             * @param[in] d Piezoelectric coupling tensor.
             *
             * @throws std::invalid_argument if `c` is not positive definite.
             * @throws std::invalid_argument if `epsilon` is not positive definite.
             * @throws std::invalid_argument if the Schur complement
             *
             * ```
             * c-dᵀϵ⁻¹d
             * ```
             *
             * is not positive definite.
             */
            LinearPiezoelectricMaterial(const StiffnessTensor &c, const PermittivityTensor &epsilon, const CouplingTensor &d)
                : c_(c), epsilon_(epsilon), d_(d) {

                if (!detail::isPD(c_)) {
                    throw std::invalid_argument("Stiffness tensor is not positive definite.");
                }

                if (!detail::isPD(epsilon_)) {
                    throw std::invalid_argument("Permittivity tensor is not positive definite.");
                }

                op_ << c_, -d_.transpose(),
                       -d_, -epsilon_;

                // Schur complement must be positive definite for thermodynamic stability.
                const auto schur = c_ - d_.transpose() * epsilon_.inverse() * d_;
                if (!detail::isPD(schur)) {
                    throw std::invalid_argument("Schur complement is not positive definite.");
                }
            }

            // /**
            //  * @brief Converts from a compatible piezoelectric material type.
            //  *
            //  * This constructor allows piezoelectric materials built from compatible
            //  * mechanical and electrical material types to be converted into this
            //  * instantiation.
            //  *
            //  * This is needed when solver code expects a canonical
            //  * `LinearPiezoelectricMaterial<...>` type, but the input uses derived
            //  * or aliased component material types.
            //  *
            //  * @tparam OtherMechanicalMaterial Compatible mechanical material type.
            //  * @tparam OtherElectricalMaterial Compatible electrical material type.
            //  *
            //  * @param[in] other Piezoelectric material to convert from.
            //  */
            // template <typename OtherMechanicalMaterial, typename OtherElectricalMaterial>
            // LinearPiezoelectricMaterial(const LinearPiezoelectricMaterial<OtherMechanicalMaterial, OtherElectricalMaterial> &other)
            //     : LinearPiezoelectricMaterial(other.elasticMaterial(), other.dielectricMaterial(), other.couplingTensor()) {}

            /// @brief Stiffness tensor.
            const StiffnessTensor &stiffnessTensor() const noexcept {
                return c_;
            }

            /// @brief Permittivity tensor.
            const PermittivityTensor &permittivityTensor() const noexcept {
                return epsilon_;
            }

            /// @brief Piezoelectric coupling tensor.
            const CouplingTensor &couplingTensor() const noexcept {
                return d_;
            }

            /** @brief Coupled constitutive operator.
             *
             * ```text
             * ⎡ c  -dᵀ⎤
             * ⎣-d  -ϵ ⎦
             * ```
             */
            const MaterialTensor &materialTensor() const noexcept {
                return op_;
            }

            /// @brief Equality comparison.
            bool operator==(const LinearPiezoelectricMaterial &other) const noexcept {
                return c_.isApprox(other.c_, NUMERICAL_ZERO)
                    && epsilon_.isApprox(other.epsilon_, NUMERICAL_ZERO)
                    && d_.isApprox(other.d_, NUMERICAL_ZERO);
            }

            /// @brief Inequality comparison.
            bool operator!=(const LinearPiezoelectricMaterial &other) const noexcept {
                return !(*this == other);
            }

        private:
            /// @brief Stiffness tensor.
            const StiffnessTensor c_;

            /// @brief Permittivity tensor.
            const PermittivityTensor epsilon_;

            /// @brief Piezoelectric coupling tensor.
            const CouplingTensor d_;

            /** @brief Coupled constitutive operator.
             *
             * ```text
             * ⎡ c  -dᵀ⎤
             * ⎣-d  -ϵ ⎦
             * ```
             */
            MaterialTensor op_;
        };

    } // namespace material

} // namespace monad
