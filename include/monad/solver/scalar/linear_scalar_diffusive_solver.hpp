#pragma once

#include "monad/solver/scalar/linear_scalar_diffusive_physics.hpp"
#include "monad/solver/periodic_cell_solver.hpp"

namespace monad {

    /**
     * @brief Periodic unit cell solver for linear scalar diffusive problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     * @tparam C Gradient sign convention (GradientConvention::Negative or Positive).
     */
    template <class Grid, class Element, fem::scalar::GradientConvention C>
    class LinearScalarDiffusiveSolver : public solver::PeriodicCellSolver<Grid, Element, solver::scalar::LinearScalarDiffusivePhysics<Element, C>> {
    public:
        using Base = solver::PeriodicCellSolver<Grid, Element, solver::scalar::LinearScalarDiffusivePhysics<Element, C>>;
        using Base::Base;
    };

} // namespace monad
