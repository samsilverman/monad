#pragma once

#include <array>
#include <cstddef>
#include <Eigen/Core>

namespace monad {

    namespace integration {

        /**
         * @brief Quadrature rule for numerical integration.
         *
         * A quadrature rule contains:
         *
         * 1. A list of integration points.
         *
         * 2. A corresponding list of weights.
         *
         * @tparam D Spatial dimension of the integration domain.
         * @tparam N Number of integration points.
         */
        template <int D, int N>
        struct QuadratureRule {
            static_assert(D > 0, "Spatial dimension of the integration domain must be positive.");
            static_assert(N > 0, "Number of integration points must be positive.");

            /// @brief Spatial dimension of the integration domain.
            static constexpr int Dim = D;

            /// @brief Number of integration points.
            static constexpr int NumPoints = N;

            using Point = Eigen::Vector<double, Dim>;
            using PointsList = std::array<Point, NumPoints>;
            using WeightsList = std::array<double, NumPoints>;

            /// @brief Integration points.
            PointsList points;

            /// @brief Integration point weights.
            WeightsList weights;

            /// @brief Equality comparison.
            bool operator==(QuadratureRule const &other) const {
                for (std::size_t i = 0; i < NumPoints; ++i) {
                    if (!points[i].isApprox(other.points[i])) {
                        return false;
                    }
                }
                return weights == other.weights;
            }

            /// @brief Inequality comparison.
            bool operator!=(QuadratureRule const &rhs) const {
                return !(*this == rhs);
            }
        };

    } // namespace integration

} // namespace monad
