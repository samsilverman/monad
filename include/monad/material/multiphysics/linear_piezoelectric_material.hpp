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
         * @tparam MechanicalMaterial Linear elastic material type (e.g. LinearElasticMaterial2d).
         * @tparam ElectricalMaterial Linear dielectric material type (e.g. LinearDielectricMaterial2d).
         */
        template <class MechanicalMaterial, class ElectricalMaterial>
        class LinearPiezoelectricMaterial {
        public:
            static_assert(MechanicalMaterial::Dim == ElectricalMaterial::Dim, "Spatial dimension of materials must be equal.");
            static_assert(MechanicalMaterial::Dim == 2 || MechanicalMaterial::Dim == 3, "Spatial dimension D must be 2 or 3.");

            /// @brief Spatial dimension (2 or 3).
            static constexpr int Dim = MechanicalMaterial::Dim;

            /// @brief Number of components in Voigt notation.
            static constexpr int VoigtSize = MechanicalMaterial::VoigtSize;

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
             * @param[in] elasticMaterial Linear elastic material.
             * @param[in] dielectricMaterial Linear dielectric material.
             * @param[in] d Piezoelectric coupling tensor.
             *
             * @throws std::invalid_argument if the Schur complement
             *
             * ```
             * c-dᵀϵ⁻¹d
             * ```
             *
             * is not positive definite.
             */
            LinearPiezoelectricMaterial(const MechanicalMaterial &elasticMaterial, const ElectricalMaterial &dielectricMaterial, const CouplingTensor &d)
                : elasticMaterial_(elasticMaterial), dielectricMaterial_(dielectricMaterial), d_(d) {
                const auto &c = elasticMaterial_.materialTensor();
                const auto &epsilon = dielectricMaterial_.materialTensor();

                // Schur complement must be positive definite for thermodynamic stability.
                const auto schur = c - d_.transpose() * epsilon.inverse() * d_;
                if (!detail::isPD(schur)) {
                    throw std::invalid_argument("Schur complement is not positive definite.");
                }

                op_ << c, -d_.transpose(),
                       -d_, -epsilon;
            }

            /**
             * @brief Converts from a compatible piezoelectric material type.
             *
             * This constructor allows piezoelectric materials built from compatible
             * mechanical and electrical material types to be converted into this
             * instantiation.
             *
             * This is needed when solver code expects a canonical
             * `LinearPiezoelectricMaterial<...>` type, but the input uses derived
             * or aliased component material types.
             *
             * @tparam OtherMechanicalMaterial Compatible mechanical material type.
             * @tparam OtherElectricalMaterial Compatible electrical material type.
             *
             * @param[in] other Piezoelectric material to convert from.
             */
            template <typename OtherMechanicalMaterial, typename OtherElectricalMaterial>
            LinearPiezoelectricMaterial(const LinearPiezoelectricMaterial<OtherMechanicalMaterial, OtherElectricalMaterial> &other)
                : LinearPiezoelectricMaterial(other.elasticMaterial(), other.dielectricMaterial(), other.couplingTensor()) {}

            /// @brief Linear elastic material.
            const MechanicalMaterial &elasticMaterial() const noexcept {
                return elasticMaterial_;
            }

            /// @brief Linear dielectric material.
            const ElectricalMaterial &dielectricMaterial() const noexcept {
                return dielectricMaterial_;
            }

            /// @brief Piezoelectric coupling tensor.
            const CouplingTensor &couplingTensor() const noexcept {
                return d_;
            }

            /// @brief Coupled constitutive operator.
            const MaterialTensor &materialTensor() const noexcept {
                return op_;
            }

            /// @brief Equality comparison.
            bool operator==(const LinearPiezoelectricMaterial &other) const noexcept {
                return op_ == other.op_;
            }

            /// @brief Inequality comparison.
            bool operator!=(const LinearPiezoelectricMaterial &other) const noexcept {
                return !(*this == other);
            }

        private:
            /// @brief Linear elastic material.
            const MechanicalMaterial elasticMaterial_;

            /// @brief Linear dielectric material.
            const ElectricalMaterial dielectricMaterial_;

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
