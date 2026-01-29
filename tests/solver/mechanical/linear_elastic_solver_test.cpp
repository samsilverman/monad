#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/solver/mechanical/linear_elastic_solver.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::detail;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::LinearElasticSolver2d: Test solve", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Material = LinearElasticMaterial<TestType::Dim>;
    using StiffnessTensor = typename Material::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);
    
    StiffnessTensor C = StiffnessTensor::Random();
    // Make PSD
    C = C.transpose() * C;
    // Make PD
    C += StiffnessTensor::Identity();

    const Material material(C);

    SECTION("Density=1 → CBar=material") {
        grid.setDensitiesOnes();
        const LinearElasticSolver solver(grid, material);

        const auto results = solver.solve();
        const auto &CBar = results.CBar;

        const auto expected = material.materialTensor();

        REQUIRE(CBar.isApprox(expected, NUMERICAL_ZERO));
    }

    SECTION("Density=0 → CBar=0") {
        grid.setDensitiesZeros();
        LinearElasticSolver solver(grid, material);

        const auto results = solver.solve();
        const auto &CBar = results.CBar;

        // Set zero threshold to one order of magnitude larger
        REQUIRE(CBar.isZero(10.0 * NUMERICAL_ZERO));
    }

    SECTION("Random density") {
        grid.setDensitiesRandom(1234);
        const LinearElasticSolver solver(grid, material);

        auto results = solver.solve();
        const auto &CBar = results.CBar;

        REQUIRE(isSymmetric(CBar));
        REQUIRE(isPD(CBar));

        SECTION("Voigt/Reuss bounds") {
            const auto voigt = material.voigt(grid);
            const auto reuss = material.reuss(grid);

            REQUIRE(reuss.trace() <= CBar.trace());
            REQUIRE(CBar.trace() <= voigt.trace());
        }

        SECTION("Translational invariance") {
            Resolution shift;
            shift.fill(1);

            grid.translate(shift);

            const LinearElasticSolver solver2(grid, material);

            results = solver.solve();
            const auto &CBar2 = results.CBar;

            REQUIRE(CBar2.isApprox(CBar, NUMERICAL_ZERO));
        }
    }

    SECTION("Bad options") {
        grid.setDensitiesRandom(1234);
        const LinearElasticSolver solver(grid, material);

        SolverOptions options;
        options.maxIterations = 1;

        REQUIRE_THROWS_AS(solver.solve(options), std::runtime_error);
    }
}
