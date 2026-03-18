#pragma once

#include "monad/field/density_field.hpp"
#include "monad/detail/mean.hpp"

namespace monad {
    
    /**
     * @brief Voigt upper bound for the homogenized material tensor.
     *
     * @tparam Material Material type (e.g. LinearElasticMaterial2d).
     *
     * @param[in] material Base material.
     * @param[in] densityField Per-element density field.
     *
     * @returns Voigt upper bound for the homogenized material tensor.
     */
    template <class Material>
    auto voigtBound(const Material& material, const field::DensityField<Material::Dim>& densityField) {
        const double density = detail::arithmeticMean(densityField.densities());

        return density * material.materialTensor();
    }

    /**
     * @brief Reuss lower bound for the homogenized material tensor.
     *
     * @tparam Material Material type (e.g. LinearElasticMaterial2d).
     *
     * @param[in] material Base material.
     * @param[in] densityField Per-element density field.
     *
     * @returns Reuss lower bound for the homogenized material tensor.
     */
    template <class Material>
    auto reussBound(const Material& material, const field::DensityField<Material::Dim>& densityField) {
        const double density = detail::harmonicMean(densityField.densities());

        return density * material.materialTensor();
    }

} // namespace monad
