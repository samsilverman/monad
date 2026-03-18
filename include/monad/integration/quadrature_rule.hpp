#pragma once

#include <array>
#include <cstddef>
#include <Eigen/Core>
#include "monad/detail/constants.hpp"

namespace monad {

    namespace integration {

        /**
         * @brief Quadrature rule for numerical integration.
         *
         * A quadrature rule consists of a set of integration
         * points and the corresponding integration weights.
         *
         * @tparam D Spatial dimension.
         * @tparam N Number of integration points.
         */
        template <int D, int N>
        struct QuadratureRule {
            static_assert(D > 0, "Spatial dimension of the integration domain must be positive.");
            static_assert(N > 0, "Number of integration points must be positive.");

            /// @brief Spatial dimension.
            static constexpr int Dim = D;

            /// @brief Number of integration points.
            static constexpr int NumPoints = N;

            using Point = Eigen::Vector<double, Dim>;
            using PointsList = std::array<Point, NumPoints>;
            using WeightsList = std::array<double, NumPoints>;

            /// @brief Integration points.
            PointsList points;

            /// @brief Integration weights.
            WeightsList weights;

            /// @brief Equality comparison.
            bool operator==(const QuadratureRule &other) const noexcept {
                for (std::size_t i = 0; i < NumPoints; ++i) {
                    if (!points[i].isApprox(other.points[i], NUMERICAL_ZERO)) {
                        return false;
                    }
                }
                return weights == other.weights;
            }

            /// @brief Inequality comparison.
            bool operator!=(const QuadratureRule &rhs) const noexcept {
                return !(*this == rhs);
            }
        };

    } // namespace integration

} // namespace monad
