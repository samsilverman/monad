#pragma once

#include <cstddef>
#include <functional>
#include <cmath>
#include <stdexcept>
#include "monad/field/density_field.hpp" 
#include "monad/integration/integrate_scalar.hpp"

namespace monad {

    /**
     * @brief Constructs a density field from a continuous function.
     *
     * Element densities are computed by averaging a continuous function
     *
     * ```text
     * f:ℝᵈ→[0,1]
     * ```
     *
     * over each grid element using numerical quadrature.
     *
     * @tparam Grid Grid type (e.g. Quad4Grid).
     *
     * @param[in] grid Grid.
     * @param[in] f Function.
     *
     * @returns Density field defined on the grid elements.
     *
     * @throws std::invalid_argument if `f` returns a value outside the range [0,1].
     *
     * @note Densities are clamped to a minimum value of `ZERO_TOLERANCE` to improve
     * numerical stability.
     */
    template <class Grid>
    inline field::DensityField<Grid::Dim> makeDensityFieldFromFunction(const Grid &grid, const std::function<double(const typename Grid::Point&)>& f) {
        field::DensityField<Grid::Dim> densityField(grid.resolution());

        for (std::size_t i = 0; i < grid.numElements(); ++i) {
            const auto nodes = grid.elementNodes(i);

            auto integrand = [&](const typename Grid::Point& point) -> double {
                const auto N = Grid::Element::shapeFunctions(point);
                const typename Grid::Point globalPoint = N.transpose() * nodes;

                const double density = f(globalPoint);
                if (density < 0.0 || density > 1.0) {
                    throw std::invalid_argument("Function value (" + std::to_string(density) + ") is outside range [0,1].");
                }

                const auto J = Grid::Element::jacobian(point, nodes);

                return density * std::abs(J.determinant());
            };

            const double density = integration::integrateScalar(integrand, Grid::Element::quadratureRule()) / Grid::Element::measure(nodes);
            densityField.setDensity(i, density);
        }

        return densityField;
    }

} // namespace monad
