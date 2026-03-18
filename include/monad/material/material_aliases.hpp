#pragma once

#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/mechanical/linear_elastic_material_3d.hpp"

namespace monad {

    /**
     * @brief 2D linear dielectric material model.
     *
     * By Gauss's law, the electric displacement D∈ℝ² is a linear function of
     * the electric field E∈ℝ²:
     *
     * ```text
     * D=ϵE=-ϵ∇φ
     * ```
     *
     * - φ∈ℝ is the electric potential.
     *
     * - ϵ∈Sym₂(ℝ) is the permittivity tensor.
     */
    using LinearDielectricMaterial2d = material::LinearTransportMaterial<2>;

    /**
     * @brief 3D linear dielectric material model.
     *
     * By Gauss's law, the electric displacement D∈ℝ³ is a linear function of
     * the electric field E∈ℝ³:
     *
     * ```text
     * D=ϵE=-ϵ∇φ
     * ```
     *
     * - φ∈ℝ is the electric potential.
     *
     * - ϵ∈Sym₃(ℝ) is the permittivity tensor.
     */
    using LinearDielectricMaterial3d = material::LinearTransportMaterial<3>;

    /**
     * @brief 2D linear electrical conductive material model.
     *
     * By Ohm's law, the current density J∈ℝ² is a linear function of
     * the electric field E∈ℝ²:
     *
     * ```text
     * J=σE=-σ∇φ
     * ```
     *
     * - φ∈ℝ is the electric potential.
     *
     * - σ∈Sym₂(ℝ) is the conductivity tensor.
     */
    using LinearElectricalConductiveMaterial2d = material::LinearTransportMaterial<2>;

    /**
     * @brief 3D linear electrical conductive material model.
     *
     * By Ohm's law, the current density J∈ℝ³ is a linear function of
     * the electric field E∈ℝ³:
     *
     * ```text
     * J=σE=-σ∇φ
     * ```
     *
     * - φ∈ℝ is the electric potential.
     *
     * - σ∈Sym₃(ℝ) is the conductivity tensor.
     */
    using LinearElectricalConductiveMaterial3d = material::LinearTransportMaterial<3>;

    /**
     * @brief 2D linear magnetic material model
     *
     * In a linear magnetic constitutive law, the magnetic flux density B∈ℝ² is a linear function of
     * the magnetic field H∈ℝ²:
     *
     * ```text
     * B=μH=-μ∇φ
     * ```
     *
     * - φ∈ℝ is the magnetic potential.
     *
     * - μ∈Sym₂(ℝ) is the permeability tensor.
     */
    using LinearMagneticMaterial2d = material::LinearTransportMaterial<2>;

    /**
     * @brief 3D linear magnetic material model
     *
     * In a linear magnetic constitutive law, the magnetic flux density B∈ℝ³ is a linear function of
     * the magnetic field H∈ℝ³:
     *
     * ```text
     * B=μH=-μ∇φ
     * ```
     *
     * - φ∈ℝ is the magnetic potential.
     *
     * - μ∈Sym₃(ℝ) is the permeability tensor.
     */
    using LinearMagneticMaterial3d = material::LinearTransportMaterial<3>;

    /**
     * @brief 2D linear mass diffusive material model.
     *
     * By Fick's law, the diffusion flux J∈ℝ² is a linear function of
     * the mass concentration gradient ∇c∈ℝ²:
     *
     * ```text
     * J=-D∇c
     * ```
     *
     * - c∈ℝ is the mass concentration.
     *
     * - D∈Sym₂(ℝ) is the diffusivity tensor.
     */
    using LinearMassDiffusiveMaterial2d = material::LinearTransportMaterial<2>;

    /**
     * @brief 3D linear mass diffusive material model.
     *
     * By Fick's law, the diffusion flux J∈ℝ³ is a linear function of
     * the mass concentration gradient ∇c∈ℝ³:
     *
     * ```text
     * J=-D∇c
     * ```
     *
     * - c∈ℝ is the mass concentration.
     *
     * - D∈Sym₃(ℝ) is the diffusivity tensor.
     */
    using LinearMassDiffusiveMaterial3d = material::LinearTransportMaterial<3>;

    /**
     * @brief 2D linear porous material model.
     *
     * By Darcy's law, the volumetric flux q∈ℝ² is a linear function of
     * the pressure gradient ∇p∈ℝ²:
     *
     * ```text
     * q=-K∇p
     * ```
     *
     * - p∈ℝ is the pressure.
     *
     * - K∈Sym₂(ℝ) is the permeability tensor.
     */
    using LinearPorousMaterial2d = material::LinearTransportMaterial<2>;

    /**
     * @brief 3D linear porous material model.
     *
     * By Darcy's law, the volumetric flux q∈ℝ³ is a linear function of
     * the pressure gradient ∇p∈ℝ³:
     *
     * ```text
     * q=-K∇p
     * ```
     *
     * - p∈ℝ is the pressure.
     *
     * - K∈Sym₃(ℝ) is the permeability tensor.
     */
    using LinearPorousMaterial3d = material::LinearTransportMaterial<3>;

    /**
     * @brief 2D linear thermal conductive material model.
     *
     * By Fourier's law, the heat flux q∈ℝ² is a linear function of
     * the temperature gradient ∇T∈ℝ²:
     *
     * ```text
     * q=-κ∇T
     * ```
     *
     * - T∈ℝ is the temperature.
     *
     * - κ∈Sym₂(ℝ) is the conductivity tensor.
     */
    using LinearThermalConductiveMaterial2d = material::LinearTransportMaterial<2>;

    /**
     * @brief 3D linear thermal conductive material model.
     *
     * By Fourier's law, the heat flux q∈ℝ³ is a linear function of
     * the temperature gradient ∇T∈ℝ³:
     *
     * ```text
     * q=-κ∇T
     * ```
     *
     * - T∈ℝ is the temperature.
     *
     * - κ∈Sym₃(ℝ) is the conductivity tensor.
     */
    using LinearThermalConductiveMaterial3d = material::LinearTransportMaterial<3>;

    /**
     * @brief 2D linear piezoelectric material model.
     *
     * In the stress-charge form, mechanical fields (stress T∈ℝ³ and strain S∈ℝ³)
     * and electrical fields (electric displacement D∈ℝ² and electric field E∈ℝ²)
     * are coupled by:
     *
     * ```text
     * T=cS-dᵀE
     * D=dT+ϵE
     * ```
     *
     * - c∈Sym₃(ℝ) is the stiffness tensor in Voigt notation.
     *
     * - ϵ∈Sym₂(ℝ) is the permittivity tensor.
     *
     * - d∈ℝ²ˣ³ is the piezoelectric coupling tensor.
     */
    using LinearPiezoelectricMaterial2d = material::LinearPiezoelectricMaterial<LinearElasticMaterial2d, LinearDielectricMaterial2d>;

    /**
     * @brief 3D linear piezoelectric material model.
     *
     * In the stress-charge form, mechanical fields (stress T∈ℝ⁶ and strain S∈ℝ⁶)
     * and electrical fields (electric displacement D∈ℝ³ and electric field E∈ℝ³)
     * are coupled by:
     *
     * ```text
     * T=cS-dᵀE
     * D=dT+ϵE
     * ```
     *
     * - c∈Sym₆(ℝ) is the stiffness tensor in Voigt notation.
     *
     * - ϵ∈Sym₃(ℝ) is the permittivity tensor.
     *
     * - d∈ℝ³ˣ⁶ is the piezoelectric coupling tensor.
     */
    using LinearPiezoelectricMaterial3d = material::LinearPiezoelectricMaterial<LinearElasticMaterial3d, LinearDielectricMaterial3d>;

} // namespace monad
