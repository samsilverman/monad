#pragma once

#include "monad/grid/grid_base.hpp"
#include "monad/solver/multiphysics/linear_piezoelectric_physics.hpp"
#include "monad/solver/periodic_cell_solver.hpp"
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/material/transport/linear_transport_material.hpp"

namespace monad {

    template <class Grid, class Element>
    class LinearPiezoelectricSolver : public solver::PeriodicCellSolver<Grid, Element, solver::multiphysics::LinearPiezoelectricPhysics<Element>> {
    public:
        using Base = solver::PeriodicCellSolver<Grid, Element, solver::multiphysics::LinearPiezoelectricPhysics<Element>>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element, typename MechanicalMaterial, typename ElectricalMaterial>
    LinearPiezoelectricSolver(const GridBase<Grid, Element> &, const LinearPiezoelectricMaterial<MechanicalMaterial, ElectricalMaterial> &) -> LinearPiezoelectricSolver<Grid, Element>;

} // namespace monad
