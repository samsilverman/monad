#include <vector>
#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "monad/detail/mean.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::detail;

TEST_CASE("monad::detail: Test arithmeticMean", "[monad]") {
    const std::vector<double> x{-1.1, 0.0, 1.1, 2.2, 3.3};

    const double expected = (-1.1 + 0.0 + 1.1 + 2.2 + 3.3) / 5.0;

    REQUIRE_THAT(arithmeticMean(x), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::detail: Test harmonicMean", "[monad]") {
    std::vector<double> x{-1.1, 1.1, 2.2, 3.3};

    SECTION("No errors") {
        const double expected = 4.0 / ((1.0 / -1.1) + (1.0 / 1.1) + (1.0 / 2.2) + (1.0 / 3.3));

        REQUIRE_THAT(harmonicMean(x), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
    }

    SECTION("value=0") {
        x.push_back(0.0);

        REQUIRE_THROWS_AS(harmonicMean(x), std::invalid_argument);
    }
}
