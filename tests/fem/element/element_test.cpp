#include <tuple>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "monad/fem/element/quad4.hpp"
#include "monad/fem/element/quad8.hpp"
#include "monad/fem/element/hex8.hpp"
#include "monad/fem/element/hex20.hpp"
#include "monad/integration/integrate_scalar.hpp"
#include "monad/detail/constants.hpp"
#include "test_common.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::integration;
using namespace monad::testing;

using Types = std::tuple<Quad4, Quad8, Hex8, Hex20>;

TEMPLATE_LIST_TEST_CASE("monad::fem::ElementBase Test localNodes", "[monad]", Types) {
    // Bounding box of [-1,1]ᵈ
    const auto nodes = TestType::localNodes();

    REQUIRE_THAT(nodes.minCoeff(), Catch::Matchers::WithinAbs(-1.0, NUMERICAL_ZERO));
    REQUIRE_THAT(nodes.maxCoeff(), Catch::Matchers::WithinAbs(1.0, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::fem::ElementBase Test shapeFunctions", "[monad]", Types) {
    SECTION("Corner nodes") {
        const auto points = TestType::localNodes();

        for (auto point : points.rowwise()) {
            const auto N = TestType::shapeFunctions(point);

            // Partition of unity
            REQUIRE_THAT(N.sum(), Catch::Matchers::WithinAbs(1.0, NUMERICAL_ZERO));

            // Kronecker-delta
            REQUIRE_THAT(N.maxCoeff(), Catch::Matchers::WithinAbs(1.0, NUMERICAL_ZERO));
        }
    }

    SECTION("Interior points") {
        const auto points = TestType::quadratureRule().points;

        for (auto point : points) {
            const auto N = TestType::shapeFunctions(point);

            // Partition of unity
            REQUIRE_THAT(N.sum(), Catch::Matchers::WithinAbs(1.0, NUMERICAL_ZERO));
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::ElementBase Test gradShapeFunctions", "[monad]", Types) {
    // Finite difference step size
    const double h = 1e-5;

    const auto points = TestType::quadratureRule().points;

    for (auto point : points) {
        auto dN = TestType::gradShapeFunctions(point);

        // Compare to numeric solution from finite (central) differences
        for (int i = 0; i < TestType::Dim; i++) {
            const auto analytic = dN.row(i).transpose();

            auto pointPlusH = point;
            pointPlusH(i) += h;

            auto pointMinusH = point;
            pointMinusH(i) -= h;

            // Use eval() to convert expression to vector
            const auto numeric = ((TestType::shapeFunctions(pointPlusH) - TestType::shapeFunctions(pointMinusH)) / (2.0 * h)).eval();

            REQUIRE(numeric.isApprox(analytic, NUMERICAL_ZERO));
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::ElementBase Test jacobian", "[monad]", Types) {
    using JacobianMatrix = typename TestType::JacobianMatrix;

    const auto points = TestType::quadratureRule().points;

    for (auto point : points) {
        // Scaled (x2) mapping → J=2I
        const auto J = TestType::jacobian(point, 2.0 * TestType::localNodes());

        const auto expected = 2.0 * JacobianMatrix::Identity();

        REQUIRE(J.isApprox(expected, NUMERICAL_ZERO));
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::ElementBase Test quadratureRule", "[monad]", Types) {
    using Point = typename TestType::Point;

    // Quadrature rule exact for degree <= pExact polynomials
    int pExact;
    if constexpr (std::is_same_v<TestType, Quad4> || std::is_same_v<TestType, Hex8>) {
        pExact = 3;
    }
    else {
        pExact = 5;
    }

    const auto rule = TestType::quadratureRule();

    for (int p = 0; p <= pExact + 1; ++p) {
        double expected;
        if constexpr (TestType::Dim == 2) {
            expected = analyticIntegral(p, p);
        }
        else {
            expected = analyticIntegral(p, p, p);
        }

        auto integrand = [p](const Point &point) -> double {
            if constexpr (TestType::Dim == 2) {
                return numericIntegrand(point, p, p);
            }
            else {
                return numericIntegrand(point, p, p, p);
            }
        };

        if (p <= pExact) {
            REQUIRE_THAT(integrateScalar(integrand, rule), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
        }
        else {
            REQUIRE_THAT(integrateScalar(integrand, rule), !Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::ElementBase Test measure", "[monad]", Types) {
    // [‑0.5,0.5]² → area = 1
    // [‑0.5,0.5]³ → volume = 1
    const double measure = TestType::measure(0.5 * TestType::localNodes());

    REQUIRE_THAT(measure, Catch::Matchers::WithinAbs(1.0, NUMERICAL_ZERO));
}

TEST_CASE("monad::fem::ElementBase Test gmshElementType", "[monad]") {
    SECTION("Quad4") {
        REQUIRE(Quad4::gmshElementType() == 3);
    }

    SECTION("Quad8") {
        REQUIRE(Quad8::gmshElementType() == 16);
    }

    SECTION("Hex8") {
        REQUIRE(Hex8::gmshElementType() == 5);
    }

    SECTION("Hex20") {
        REQUIRE(Hex20::gmshElementType() == 17);
    }
}

TEST_CASE("monad::fem::ElementBase Test gmshNodeOrdering", "[monad]") {
    SECTION("Quad4") {
        const auto expected = Quad4::NodeIndicesList{0, 1, 2, 3};

        REQUIRE(Quad4::gmshNodeOrdering() == expected);
    }

    SECTION("Quad8") {
        const auto expected = Quad8::NodeIndicesList{0, 1, 2, 3, 4, 5, 6, 7};

        REQUIRE(Quad8::gmshNodeOrdering() == expected);
    }

    SECTION("Hex8") {
        const auto expected = Hex8::NodeIndicesList{0, 1, 5, 4, 3, 2, 6, 7};

        REQUIRE(Hex8::gmshNodeOrdering() == expected);
    }

    SECTION("Hex20") {
        const auto expected = Hex20::NodeIndicesList{0, 1, 5, 4, 3, 2, 6, 7, 8, 16, 11, 17, 9, 12, 13, 15, 10, 19, 18, 14};

        REQUIRE(Hex20::gmshNodeOrdering() == expected);
    }
}
