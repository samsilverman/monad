#include <catch2/catch_test_macros.hpp>
#include <Eigen/Core>
#include "monad/integration/quadrature_rule.hpp"

using namespace monad::integration;

TEST_CASE("monad::integration::QuadratureRule: Test operator==", "[monad]") {
    QuadratureRule<2, 3> rule1;
    rule1.points = {
        Eigen::Vector2d(1.0, 0.4),
        Eigen::Vector2d(0.3, 0.0),
        Eigen::Vector2d(0.1, 0.1)
    };
    rule1.weights = {-0.1, 0.7, 1.1};
    
    QuadratureRule<2, 3> rule2;
    rule2.points = rule1.points;
    rule2.weights = rule1.weights;

    REQUIRE(rule1 == rule2);
}

TEST_CASE("monad::integration::QuadratureRule: Test operator!=", "[monad]") {
    QuadratureRule<2, 3> rule1;
    rule1.points = {
        Eigen::Vector2d(1.0, 0.4),
        Eigen::Vector2d(0.3, 0.0),
        Eigen::Vector2d(0.1, 0.1)
    };
    rule1.weights = {-0.1, 0.7, 1.1};

    QuadratureRule<2, 3> rule2;

    SECTION("Different points") {
        rule2.points = {
            Eigen::Vector2d(1.0, 0.4),
            Eigen::Vector2d(-0.3, 0.0),
            Eigen::Vector2d(0.1, 0.1)
        };
        rule2.weights = rule1.weights;

        REQUIRE(rule1 != rule2);
    }

    SECTION("Different weights") {
        rule2.points = rule2.points;
        rule1.weights = {-0.1, -0.7, 1.1};

        REQUIRE(rule1 != rule2);
    }
}
