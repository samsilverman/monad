#pragma once

#include "monad/grid/grid_base.hpp"
#include "monad/solver/scalar/linear_scalar_diffusive_solver.hpp"
#include "monad/material/transport/linear_transport_material.hpp"

namespace monad {

    /**
     * @brief Periodic unit cell solver for linear dielectric problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearDielectricSolver : public LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Negative> {
    public:
        using Base = LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Negative>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearDielectricSolver(const GridBase<Grid, Element> &, const LinearTransportMaterial<Element::Dim> &) -> LinearDielectricSolver<Grid, Element>;

    /**
     * @brief Periodic unit cell solver for linear electrical conductive problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearElectricalConductiveSolver : public LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Negative> {
    public:
        using Base = LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Negative>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearElectricalConductiveSolver(const GridBase<Grid, Element> &, const LinearTransportMaterial<Element::Dim> &) -> LinearElectricalConductiveSolver<Grid, Element>;

    /**
     * @brief Periodic unit cell solver for linear magnetic problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearMagneticSolver : public LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Negative> {
    public:
        using Base = LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Negative>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearMagneticSolver(const GridBase<Grid, Element> &, const LinearTransportMaterial<Element::Dim> &) -> LinearMagneticSolver<Grid, Element>;

    /**
     * @brief Periodic unit cell solver for linear mass diffusive problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearMassDiffusiveSolver : public LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Positive> {
    public:
        using Base = LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Positive>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearMassDiffusiveSolver(const GridBase<Grid, Element> &, const LinearTransportMaterial<Element::Dim> &) -> LinearMassDiffusiveSolver<Grid, Element>;

    /**
     * @brief Periodic unit cell solver for linear porous problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearPorousSolver : public LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Positive> {
    public:
        using Base = LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Positive>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearPorousSolver(const GridBase<Grid, Element> &, const LinearTransportMaterial<Element::Dim> &) -> LinearPorousSolver<Grid, Element>;

    /**
     * @brief Periodic unit cell solver for linear thermal conductive problems.
     *
     * @tparam Grid Grid class (e.g. Quad4Grid).
     * @tparam Element Element class (e.g. Quad4).
     */
    template <class Grid, class Element>
    class LinearThermalConductiveSolver : public LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Positive> {
    public:
        using Base = LinearScalarDiffusiveSolver<Grid, Element, fem::scalar::GradientConvention::Positive>;
        using Base::Base;
    };

    // Template arugment deduction
    template <class Grid, class Element>
    LinearThermalConductiveSolver(const GridBase<Grid, Element> &, const LinearTransportMaterial<Element::Dim> &) -> LinearThermalConductiveSolver<Grid, Element>;

} // namespace monad
