#pragma once

#include <type_traits>
#include <cstddef>
#include "monad/integration/quadrature_rule.hpp"

namespace monad {

    namespace integration {

        /**
         * @brief Performs numerical integration (quadrature) of a scalar integrand function.
         *
         * The integrand is a scalar function from ℝᵈ→ℝ.
         *
         * @tparam F Type of the integrand.
         * @tparam D Spatial dimension of the integration domain.
         * @tparam N Number of integration points.
         *
         * @param[in] integrand Function to be integrated.
         * @param[in] rule Quadrature rule.
         *
         * @returns Approximate value of the integral.
         */
        template <typename F, int D, int N>
        double integrateScalar(const F &integrand, const QuadratureRule<D, N> &rule) {
            using DomainType = typename QuadratureRule<D, N>::Point;

            static_assert(std::is_invocable_r_v<double, F, DomainType>, "Integrand must be callable as double(Point).");

            double result = 0.0;

            for (std::size_t i = 0; i < N; ++i) {
                const double weight = rule.weights[i];
                const auto &point = rule.points[i];

                result += weight * integrand(point);
            }

            return result;
        }

    } // namespace integration

} // namespace monad
