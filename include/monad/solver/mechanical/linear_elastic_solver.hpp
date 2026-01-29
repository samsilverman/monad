#pragma once

#include "monad/grid/grid_base.hpp"
#include "monad/solver/mechanical/linear_elastic_physics.hpp"
#include "monad/solver/periodic_cell_solver.hpp"
#include "monad/material/mechanical/linear_elastic_material.hpp"

namespace monad {

    /**
     * @brief Periodic unit cell solver for linear elastic problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearElasticSolver : public solver::PeriodicCellSolver<Grid, Element, solver::mechanical::LinearElasticPhysics<Element>> {
    public:
        using Base = solver::PeriodicCellSolver<Grid, Element, solver::mechanical::LinearElasticPhysics<Element>>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearElasticSolver(const GridBase<Grid, Element> &, const LinearElasticMaterial<Element::Dim> &) -> LinearElasticSolver<Grid, Element>;

} // namespace monad
