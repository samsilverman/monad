#pragma once

#include "monad/material/transport/linear_transport_material.hpp"

namespace monad {

    /**
     * @brief Represents a 2D linear dielectric material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For dielectric materials, Gauss's law provides the constitutive relation:
     *
     * D=εE
     *
     * - D∈ℝ² is the electric displacement (analogous to flux J)..
     *
     * - E=-∇φ∈ℝ² is the electric field (analogous to scalar potential gradient ∇φ)..
     *
     * - φ∈ℝ is the electric potential (analogous to scalar potential φ)..
     *
     * - ε∈Sym₂(ℝ) is the permittivity tensor (analogous to transport tensor K)..
     */
    using LinearDielectricMaterial2d = LinearTransportMaterial<2>;

    /**
     * @brief Represents a 3D linear dielectric material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For dielectric materials, Gauss's law provides the constitutive relation:
     *
     * D=εE
     *
     * - D∈ℝ³ is the electric displacement (analogous to flux J).
     *
     * - E=-∇φ∈ℝ³ is the electric field (analogous to scalar potential gradient ∇φ).
     *
     * - φ∈ℝ is the electric potential (analogous to scalar potential φ).
     *
     * - ε∈Sym₃(ℝ) is the permittivity tensor (analogous to transport tensor K).
     */
    using LinearDielectricMaterial3d = LinearTransportMaterial<3>;

    /**
     * @brief Represents a 2D linear electrical conductive material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For electrical conductive materials, Ohm's law provides the constitutive relation:
     *
     * J=σE
     *
     * - J∈ℝ² is the current density (analogous to flux J).
     *
     * - E=-∇φ∈ℝ² is the electric field (analogous to scalar potential gradient ∇φ).
     *
     * - φ∈ℝ is the electric potential (analogous to scalar potential φ).
     *
     * - σ∈Sym₂(ℝ) is the conductivity tensor (analogous to transport tensor K).
     */
    using LinearElectricalConductiveMaterial2d = LinearTransportMaterial<2>;

    /**
     * @brief Represents a 3D linear electrical conductive material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For electrical conductive materials, Ohm's law provides the constitutive relation:
     *
     * J=σE
     *
     * - J∈ℝ³ is the current density (analogous to flux J).
     *
     * - E=-∇φ∈ℝ³ is the electric field (analogous to scalar potential gradient ∇φ).
     *
     * - φ∈ℝ is the electric potential (analogous to scalar potential φ).
     *
     * - σ∈Sym₃(ℝ) is the conductivity tensor (analogous to transport tensor K).
     */
    using LinearElectricalConductiveMaterial3d = LinearTransportMaterial<3>;

    /**
     * @brief Represents a 2D linear magnetic material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For magnetic materials, the constitutive relation is:
     *
     * B=μH
     *
     * - B∈ℝ² is the magnetic flux density (analogous to flux J).
     *
     * - H=-∇φ∈ℝ² is the magnetic field (analogous to scalar potential gradient ∇φ).
     *
     * - φ∈ℝ is the magnetic potential (analogous to scalar potential φ).
     *
     * - μ∈Sym₂(ℝ) is the permeability tensor (analogous to transport tensor K).
     */
    using LinearMagneticMaterial2d = LinearTransportMaterial<2>;

    /**
     * @brief Represents a 3D linear magnetic material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For magnetic materials, the constitutive relation is:
     *
     * B=μH
     *
     * - B∈ℝ³ is the magnetic flux density (analogous to flux J).
     *
     * - H=-∇φ∈ℝ³ is the magnetic field (analogous to scalar potential gradient ∇φ).
     *
     * - φ∈ℝ is the magnetic potential (analogous to scalar potential φ).
     *
     * - μ∈Sym₃(ℝ) is the permeability tensor (analogous to transport tensor K).
     */
    using LinearMagneticMaterial3d = LinearTransportMaterial<3>;

    /**
     * @brief Represents a 2D linear mass diffusive material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For mass diffusive materials, Fick's law provides the constitutive relation:
     *
     * J=-D∇c
     *
     * - J∈ℝ² is the diffusion flux (analogous to flux J).
     *
     * - ∇c∈ℝ² is the mass concentration gradient (analogous to scalar potential gradient ∇φ).
     *
     * - D∈Sym₂(ℝ) is the diffusivity tensor (analogous to transport tensor K).
     */
    using LinearMassDiffusiveMaterial2d = LinearTransportMaterial<2>;

    /**
     * @brief Represents a 3D linear mass diffusive material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For mass diffusive materials, Fick's law provides the constitutive relation:
     *
     * J=-D∇c
     *
     * - J∈ℝ³ is the diffusion flux (analogous to flux J).
     *
     * - ∇c∈ℝ³ is the mass concentration gradient (analogous to scalar potential gradient ∇φ).
     *
     * - D∈Sym₃(ℝ) is the diffusivity tensor (analogous to transport tensor K).
     */
    using LinearMassDiffusiveMaterial3d = LinearTransportMaterial<3>;

    /**
     * @brief Represents a 2D linear porous material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For porous materials, Darcy's law provides the constitutive relation:
     *
     * q=-K∇p
     *
     * - q∈ℝ² is the volumetric flux (analogous to flux J).
     *
     * - ∇p∈ℝ² is the pressure gradient (analogous to scalar potential gradient ∇φ).
     *
     * - K∈Sym₂(ℝ) is the permeability tensor (analogous to transport tensor K).
     */
    using LinearPorousMaterial2d = LinearTransportMaterial<2>;

    /**
     * @brief Represents a 3D linear porous material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For porous materials, Darcy's law provides the constitutive relation:
     *
     * q=-K∇p
     *
     * - D∈ℝ³ is the volumetric flux (analogous to flux J).
     *
     * - ∇p∈ℝ³ is the pressure gradient (analogous to scalar potential gradient ∇φ).
     *
     * - K∈Sym₃(ℝ) is the permeability tensor (analogous to transport tensor K).
     */
    using LinearPorousMaterial3d = LinearTransportMaterial<3>;

    /**
     * @brief Represents a 2D linear thermal conductive material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For thermal conductive materials, Fourier's law provides the constitutive relation:
     *
     * q=-κ∇T
     *
     * - q∈ℝ² is the heat flux (analogous to flux J).
     *
     * - ∇T∈ℝ² is the temperature gradient (analogous to scalar potential gradient ∇φ).
     *
     * - κ∈Sym₂(ℝ) is the conductivity tensor (analogous to transport tensor K).
     */
    using LinearThermalConductiveMaterial2d = LinearTransportMaterial<2>;

    /**
     * @brief Represents a 3D linear thermal conductive material model.
     *
     * This type is a physics-specific alias of a general linear transport material,
     * whose constitutive relation takes the form:
     *
     * J=-K∇φ
     *
     * For thermal conductive materials, Fourier's law provides the constitutive relation:
     *
     * q=-κ∇T
     *
     * - q∈ℝ³ is the heat flux (analogous to flux J).
     *
     * - ∇T∈ℝ³ is the temperature gradient (analogous to scalar potential gradient ∇φ).
     *
     * - κ∈Sym₃(ℝ) is the conductivity tensor (analogous to transport tensor K).
     */
    using LinearThermalConductiveMaterial3d = LinearTransportMaterial<3>;

} // namespace monad
