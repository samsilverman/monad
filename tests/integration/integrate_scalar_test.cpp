#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Eigen/Core>
#include "monad/integration/quadrature_rule.hpp"
#include "monad/integration/integrate_scalar.hpp"
#include "monad/detail/constants.hpp"
#include "test_common.hpp"

using namespace monad;
using namespace monad::integration;
using namespace monad::testing;

TEST_CASE("monad::integration::integrateScalar: Test integrateScalar", "[monad]") {
    // 1-point 2D Guassian quadrature rule
    QuadratureRule<2, 1> rule;
    rule.points = {
        Eigen::Vector2d::Zero()
    };
    rule.weights = { 4 };

    // 1-point rule exact for degree <= 1 polynomials
    // Use ∬xᵃyᵃdxdy for x,y∈[-1,1]
    for (int a = 0; a <= 2; ++a) {
        const double expected = analyticIntegral(a, a);

        auto integrand = [&a](const Eigen::Vector2d &point) -> double {
            return numericIntegrand(point, a, a);
        };

        if (a < 2) {
            REQUIRE_THAT(integrateScalar(integrand, rule), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
        }
        else {
            REQUIRE_THAT(integrateScalar(integrand, rule), !Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
        }
    }
}
