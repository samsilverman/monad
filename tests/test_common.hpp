#pragma once

#include <cmath>
#include <Eigen/Core>

namespace monad {

    namespace testing {

        /// @brief The analytical integral to ∬xᵃyᵇdxdy for x,y∈[-1,1].
        inline double analyticIntegral(int a, int b) noexcept {
            if ((a % 2) || (b % 2)) {
                return 0;
            }

            return 4.0 / ((a + 1) * (b + 1));
        }

        /// @brief The analytical integral to ∭xᵃyᵇzᶜdxdydz for x,y,z∈[-1,1].
        inline double analyticIntegral(int a, int b, int c) noexcept {
            if ((a % 2) || (b % 2) || (c % 2)) {
                return 0;
            }

            return 8.0 / ((a + 1) * (b + 1) * (c + 1));
        }

        /// @brief The integrand f(x,y)=xᵃyᵇ.
        inline double numericIntegrand(const Eigen::Vector2d &point, int a, int b) noexcept {
            const double x = point(0);
            const double y = point(1);

            return std::pow(x, a) * std::pow(y, b);
        }

        /// @brief The integrand f(x,y,z)=xᵃyᵇzᶜ.
        inline double numericIntegrand(const Eigen::Vector3d &point, int a, int b, int c) noexcept {
            const double x = point(0);
            const double y = point(1);
            const double z = point(2);

            return std::pow(x, a) * std::pow(y, b) * std::pow(z, c);
        }

    } // namespace testing

} // namespace monad
