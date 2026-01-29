#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/solver/scalar/linear_scalar_diffusive_solver.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::detail;
using namespace monad::fem::scalar;

template <class GridT, GradientConvention C>
struct TypePair {
    using Grid = GridT;
    static constexpr GradientConvention Convention = C;
};

using Types = std::tuple<
    TypePair<Quad4Grid, GradientConvention::Negative>,
    TypePair<Quad4Grid, GradientConvention::Positive>,
    TypePair<Quad8Grid, GradientConvention::Negative>,
    TypePair<Quad8Grid, GradientConvention::Positive>,
    TypePair<Hex8Grid, GradientConvention::Negative>,
    TypePair<Hex8Grid, GradientConvention::Positive>,
    TypePair<Hex20Grid, GradientConvention::Negative>,
    TypePair<Hex20Grid, GradientConvention::Positive>
>;

TEMPLATE_LIST_TEST_CASE("monad::Solver: Test solve", "[monad]", Types) {
    using Grid = typename TestType::Grid;
    using Resolution = typename Grid::Resolution;
    using Size = typename Grid::Size;
    using Element = typename Grid::Element;
    using Solver = LinearScalarDiffusiveSolver<Grid, Element, TestType::Convention>;

    Resolution resolution;
    resolution.fill(3);
    
    Size size;
    size.fill(0.5);

    Grid grid(resolution, size);

    const LinearTransportMaterial<Element::Dim> material(2.1);

    SECTION("Density=1 → KBar=material") {
        grid.setDensitiesOnes();
        const Solver solver(grid, material);

        const auto results = solver.solve();
        const auto &KBar = results.KBar;

        const auto expected = material.materialTensor();

        REQUIRE(KBar.isApprox(expected, NUMERICAL_ZERO));
    }

    SECTION("Density=0 → KBar=0") {
        grid.setDensitiesZeros();
        const Solver solver(grid, material);

        const auto results = solver.solve();
        const auto &KBar = results.KBar;

        // Set zero threshold to one order of magnitude larger
        REQUIRE(KBar.isZero(10.0 * NUMERICAL_ZERO));
    }

    SECTION("Random density") {
        grid.setDensitiesRandom(1234);
        const Solver solver(grid, material);

        auto results = solver.solve();
        const auto &KBar = results.KBar;

        REQUIRE(isSymmetric(KBar));
        REQUIRE(isPD(KBar));

        SECTION("Voigt/Reuss bounds") {
            const auto voigt = material.voigt(grid);
            const auto reuss = material.reuss(grid);

            REQUIRE(reuss.trace() <= KBar.trace());
            REQUIRE(KBar.trace() <= voigt.trace());
        }

        SECTION("Translational invariance") {
            Resolution shift;
            shift.fill(1);

            grid.translate(shift);

            const Solver solver2(grid, material);

            results = solver.solve();
            const auto &KBar2 = results.KBar;

            REQUIRE(KBar2.isApprox(KBar, NUMERICAL_ZERO));
        }
    }

    SECTION("Bad options") {
        grid.setDensitiesRandom(1234);
        const Solver solver(grid, material);

        SolverOptions options;
        options.maxIterations = 1;

        REQUIRE_THROWS_AS(solver.solve(options), std::runtime_error);
    }
}
