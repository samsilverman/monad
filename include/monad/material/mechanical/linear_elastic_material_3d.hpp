#pragma once

#include "monad/material/mechanical/linear_elastic_material.hpp"

namespace monad {

    /**
     * @brief Represents a 3D linear elastic material model.
     *
     * By Hooke's law, stress σ∈ℝ³ is a linear function of strain ε∈ℝ³:
     *
     * σ=Cε
     *
     * - C∈Sym₃(ℝ) is the stiffness tensor (in Voigt notation).
     */
    class LinearElasticMaterial3d : public LinearElasticMaterial<3> {
    public:
        using LinearElasticMaterial<3>::LinearElasticMaterial;

        /**
         * @brief Constructs a 3D isotropic linear elastic material.
         *
         * @param[in] E Young's modulus.
         * @param[in] nu Poisson's ratio.
         *
         * @throws std::invalid_argument if `E` is non-positive.
         * @throws std::invalid_argument if `nu` is not in range (-1,0.5).
         */
        LinearElasticMaterial3d(double E, double nu);
    };

} // namespace monad
