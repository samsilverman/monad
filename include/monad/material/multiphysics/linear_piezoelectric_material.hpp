#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/grid/grid_base.hpp"

namespace monad {

    /**
     * @brief Represents a linear piezoelectric material model.
     *
     * By the stress-charge form, mechanical (stress S∈ℝᵛ and strain T∈ℝᵛ) and
     * electrical (electric displacement D∈ℝᵈ and electric field E∈ℝᵈ) fields
     * are coupled:
     *
     * S=cT-dᵀE
     *
     * -D=-dT-εE
     *
     * - c∈Symᵥ(ℝ) is the stiffness tensor (in Voigt notation).
     *
     * - ε∈Symᴅ(ℝ) is the permittivity tensor (analogous to transport tensor K).
     *
     * - d∈ℝᵈˣᵛ is the piezoelectric coupling tensor.
     *
     * @tparam MechanicalMaterial Linear elastic material class (e.g. LinearElasticMaterial2d).
     * @tparam ElectricalMaterial Linear dielectric material class (e.g. LinearDielectricMaterial2d).
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
         * ⎣-d  -ε ⎦
         * ```
         */
        using MaterialTensor = Eigen::Matrix<double, VoigtSize + Dim, VoigtSize + Dim>;

        /**
         * @brief Constructs a linear piezoelectric material.
         *
         * @param[in] elasticMaterial The linear elastic material.
         * @param[in] dielectricMaterial The linear dielectric material.
         * @param[in] d The piezoelectric coupling tensor.
         *
         * @throws std::invalid_argument if Schur complement
         *
         * c-dᵀε⁻¹d
         *
         * is not positive definite.
         */
        LinearPiezoelectricMaterial(const MechanicalMaterial &elasticMaterial, const ElectricalMaterial &dielectricMaterial, const CouplingTensor &d)
            : elasticMaterial_(elasticMaterial), dielectricMaterial_(dielectricMaterial), d_(d) {
            const auto &c = elasticMaterial_.materialTensor();
            const auto &epsilon = dielectricMaterial_.materialTensor();

            // Schur complement
            // Must be PD for thermodynamic stability
            auto S = c - d_.transpose() * epsilon.inverse() * d_;
            if (!detail::isPD(S)) {
                throw std::invalid_argument("Schur complement is not positive definite.");
            }

            op_ << c, -d_.transpose(),
                   -d_, -epsilon;
        }

        /**
         * @brief Converting constructor for compatible piezoelectric material types.
         *
         * Converts any derived instantiation of `LinearPiezoelectricMaterial` into the
         * canonical form expected by the solver.
         *
         * Example:
         *
         * LinearPiezoelectricMaterial2d<LinearElasticMaterial2d, LinearDielectricMaterial2d>
         *
         * ↳ LinearPiezoelectricMaterial<LinearElasticMaterial<2>, LinearTransportMaterial<2>>
         *
         * @tparam OtherMechanicalMaterial Derived mechanical material type.
         * @tparam OtherElectricalMaterial Derived electrical material type.
         *
         * @param[in] other Linear piezoelectric material.
         *
         * @note Unlike single-physics materials, the piezoelectric class is a nested
         * template over two material types, which prevents implicit type deduction
         * for the solver. This constructor explicitly bridges that gap by forwarding
         * the material to the canonical base type.
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

        /// @brief Piezoelectric coupling tensor d.
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
         * ⎣-d  -ε ⎦
         * ```
         */
        MaterialTensor op_;
    };

    /**
     * @brief Represents a 2D linear piezoelectric material model.
     *
     * By the stress-charge form, mechanical (stress S∈ℝ³ and strain T∈ℝ³) and
     * electrical (electric displacement D∈ℝ² and electric field E∈ℝ²) fields
     * are coupled:
     *
     * S=cT-dᵀE
     *
     * -D=-dT-εE
     *
     * - c∈Sym₃(ℝ) is the stiffness tensor (in Voigt notation).
     *
     * - ε∈Sym₂(ℝ) is the dielectric tensor.
     *
     * - d∈ℝ²ˣ³ is the piezoelectric coupling tensor.
     */
    using LinearPiezoelectricMaterial2d = LinearPiezoelectricMaterial<LinearElasticMaterial<2>, LinearTransportMaterial<2>>;

    /**
     * @brief Represents a 3D linear piezoelectric material model.
     *
     * By the stress-charge form, mechanical (stress S∈ℝ⁶ and strain T∈ℝ⁶) and
     * electrical (electric displacement D∈ℝ³ and electric field E∈ℝ³) fields
     * are coupled:
     *
     * S=cT-dᵀE
     *
     * -D=-dT-εE
     *
     * - c∈Sym₆(ℝ) is the stiffness tensor (in Voigt notation).
     *
     * - ε∈Sym₃(ℝ) is the dielectric tensor.
     *
     * - d∈ℝ³ˣ⁶ is the piezoelectric coupling tensor.
     */
    using LinearPiezoelectricMaterial3d = LinearPiezoelectricMaterial<LinearElasticMaterial<3>, LinearTransportMaterial<3>>;

} // namespace monad
