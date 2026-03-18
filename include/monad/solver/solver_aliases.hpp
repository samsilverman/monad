#pragma once

#include "monad/solver/mechanical/linear_elastic_policy.hpp"
#include "monad/solver/scalar/linear_scalar_diffusive_policy.hpp"
#include "monad/solver/multiphysics/linear_piezoelectric_policy.hpp"
#include "monad/solver/periodic_cell_solver.hpp"

namespace monad {
    
    /**
     * @brief Stateless periodic-cell solver for structured linear elastic homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearElasticSolver = solver::PeriodicCellSolver<Grid, solver::mechanical::LinearElasticPolicy<typename Grid::Element>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear dielectric homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearDielectricSolver = solver::PeriodicCellSolver<Grid, solver::scalar::LinearScalarDiffusivePolicy<typename Grid::Element, fem::scalar::GradientConvention::Negative>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear electrical conductive homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearElectricalConductiveSolver = solver::PeriodicCellSolver<Grid, solver::scalar::LinearScalarDiffusivePolicy<typename Grid::Element, fem::scalar::GradientConvention::Negative>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear magnetic homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearMagneticSolver = solver::PeriodicCellSolver<Grid, solver::scalar::LinearScalarDiffusivePolicy<typename Grid::Element, fem::scalar::GradientConvention::Negative>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear mass diffusive homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearMassDiffusiveSolver = solver::PeriodicCellSolver<Grid, solver::scalar::LinearScalarDiffusivePolicy<typename Grid::Element, fem::scalar::GradientConvention::Positive>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear porous homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearPorousSolver = solver::PeriodicCellSolver<Grid, solver::scalar::LinearScalarDiffusivePolicy<typename Grid::Element, fem::scalar::GradientConvention::Negative>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear thermal conductive homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearThermalConductiveSolver = solver::PeriodicCellSolver<Grid, solver::scalar::LinearScalarDiffusivePolicy<typename Grid::Element, fem::scalar::GradientConvention::Positive>>;

    /**
     * @brief Stateless periodic-cell solver for structured linear piezoelectric homogenization problems.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     */
    template <class Grid>
    using LinearPiezoelectricSolver = solver::PeriodicCellSolver<Grid, solver::multiphysics::LinearPiezoelectricPolicy<typename Grid::Element>>;

} // namespace monad
