#pragma once

#include "monad/material/mechanical/linear_elastic_material.hpp"

namespace monad {

    /**
     * @brief 2D linear elastic material model.
     *
     * By Hooke's law, the stress σ∈ℝ³ is a linear function of
     * the strain ε∈ℝ³:
     *
     * ```text
     * σ=Cε
     * ```
     *
     * - C∈Sym₃(ℝ) is the stiffness tensor in Voigt notation.
     */
    class LinearElasticMaterial2d : public material::LinearElasticMaterial<2> {
    public:
        using material::LinearElasticMaterial<2>::LinearElasticMaterial;

        /// @brief Plane condition for 2D linear elasticity.
        enum class PlaneCondition {
            /// @brief Zero out-of-plane stress.
            PlaneStress = 0,

            /// @brief Zero out-of-plane strain.
            PlaneStrain
        };

        /**
         * @brief Constructs a 2D isotropic linear elastic material.
         *
         * @param[in] E Young's modulus.
         * @param[in] nu Poisson's ratio.
         * @param[in] condition Plane condition.
         *
         * @throws std::invalid_argument if `E` is non-positive.
         * @throws std::invalid_argument if `nu` is not in range (-1,0.5).
         */
        LinearElasticMaterial2d(double E, double nu, PlaneCondition condition);
    };

} // namespace monad
