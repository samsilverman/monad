#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/field/make_density_field_from_function.hpp"

using namespace monad;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad: Test makeDensityFieldFromFunction", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Point = typename TestType::Point;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    SECTION("No errors") {
        // f(x)=0.1Σᵢxᵢ
        auto f = [](const Point &x) -> double {
            return 0.1 * x.sum();
        };

        const auto densityField = makeDensityFieldFromFunction(grid, f);

        // The average of a linear function over a convex region equals its value at the centroid
        const auto nodes = grid.elementNodes(1);
        const Point centroid = nodes.colwise().mean();
        const double expected = f(centroid);

        REQUIRE_THAT(densityField.getDensity(1), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
    }

    SECTION("Bad function") {
        // f(x)=exp(Σᵢxᵢ)
        auto f = [](const Point &x) -> double {
            return std::exp(x.sum());
        };

        REQUIRE_THROWS_AS(makeDensityFieldFromFunction(grid, f), std::invalid_argument);
    }
}
