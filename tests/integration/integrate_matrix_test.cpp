#include <catch2/catch_test_macros.hpp>
#include <Eigen/Core>
#include "monad/integration/quadrature_rule.hpp"
#include "monad/integration/integrate_matrix.hpp"
#include "monad/detail/constants.hpp"
#include "test_common.hpp"

using namespace monad;
using namespace monad::integration;
using namespace monad::testing;

TEST_CASE("monad::integration::integrateMatrix: Test integrateMatrix", "[monad]") {
    // 1-point 2D Guassian quadrature rule
    QuadratureRule<2, 1> rule;
    rule.points = {
        Eigen::Vector2d::Zero()
    };
    rule.weights = {4};

    // 1-point rule exact for degree <= 1 polynomials
    // Use ∬xᵃyᵃdxdy for x,y∈[-1,1]
    for (int a = 0; a <= 2; ++a) {
       Eigen::Vector2d analytic = Eigen::Vector2d::Constant(analyticIntegral(a, a));

        auto integrand = [&a](const Eigen::Vector2d &point) -> Eigen::Vector2d {
            return Eigen::Vector2d::Constant(numericIntegrand(point, a, a));
        };

        const auto numeric = integrateMatrix(integrand, rule);

        if (a < 2) {
            REQUIRE(numeric.isApprox(analytic, NUMERICAL_ZERO));
        }
        else {
            REQUIRE(!numeric.isApprox(analytic, NUMERICAL_ZERO));
        }
    }
}
