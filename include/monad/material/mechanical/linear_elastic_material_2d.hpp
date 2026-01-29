#pragma once

#include "monad/material/mechanical/linear_elastic_material.hpp"

namespace monad {

    /**
     * @brief Represents a 2D linear elastic material model.
     *
     * By Hooke's law, stress σ∈ℝ² is a linear function of strain ε∈ℝ²:
     *
     * σ=Cε
     *
     * - C∈Sym₂(ℝ) is the stiffness tensor (in Voigt notation).
     */
    class LinearElasticMaterial2d : public LinearElasticMaterial<2> {
    public:
        using LinearElasticMaterial<2>::LinearElasticMaterial;

        /// @brief Defines the plane conditions for 2D linear elasticity.
        enum class PlaneCondition {
            /// @brief Assumes zero out-of-plane stress.
            PlaneStress = 0,

            /// @brief Assumes zero out-of-plane strain.
            PlaneStrain
        };

        /**
         * @brief Constructs a 2D isotropic linear elastic material.
         *
         * @param[in] E Young's modulus.
         * @param[in] nu Poisson's ratio.
         * @param[in] condition 2D plane condition.
         *
         * @throws std::invalid_argument if `E` is non-positive.
         * @throws std::invalid_argument if `nu` is not in range (-1,0.5).
         */
        LinearElasticMaterial2d(double E, double nu, PlaneCondition condition);
    };

} // namespace monad
