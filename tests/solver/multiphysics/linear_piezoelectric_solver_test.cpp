#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/solver/multiphysics/linear_piezoelectric_solver.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::detail;

// using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;
using Types = std::tuple<Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricSolver: Test solve", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Element = typename TestType::Element;
    using MechanicalMaterial = LinearElasticMaterial<Element::Dim>;
    using ElectricalMaterial = LinearTransportMaterial<Element::Dim>;
    using Material = LinearPiezoelectricMaterial<MechanicalMaterial, ElectricalMaterial>;
    using StiffnessTensor = typename MechanicalMaterial::MaterialTensor;
    using CouplingTensor = typename Material::CouplingTensor;
    using MaterialTensor = typename Material::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);

    StiffnessTensor c = StiffnessTensor::Random();
    // Make PSD
    c = c.transpose() * c;
    // Make PD
    c += StiffnessTensor::Identity();

    const MechanicalMaterial elasticMaterial(c);
    const ElectricalMaterial dielectricMaterial(2.1);

    const CouplingTensor d = 0.1 * CouplingTensor::Random();

    const Material material(elasticMaterial, dielectricMaterial, d);

    SECTION("Density=1 → homogenized=material") {
        grid.setDensitiesOnes();
        const LinearPiezoelectricSolver solver(grid, material);

        const auto results = solver.solve();
        const auto &cBar = results.cBar;
        const auto &epsilonBar = results.epsilonBar;
        const auto &dBar = results.dBar;

        MaterialTensor op;
        op << cBar, -dBar.transpose(),
              -dBar, -epsilonBar;

        const MaterialTensor expected = material.materialTensor();

        REQUIRE(op.isApprox(expected, NUMERICAL_ZERO));
    }

    SECTION("Density=0 → operator=0") {
        grid.setDensitiesZeros();
        const LinearPiezoelectricSolver solver(grid, material);

        const auto results = solver.solve();
        const auto &cBar = results.cBar;
        const auto &epsilonBar = results.epsilonBar;
        const auto &dBar = results.dBar;

        // Set zero threshold to one order of magnitude larger
        REQUIRE(cBar.isZero(10.0 * NUMERICAL_ZERO));
        REQUIRE(epsilonBar.isZero(10.0 * NUMERICAL_ZERO));
        REQUIRE(dBar.isZero(10.0 * NUMERICAL_ZERO));
    }

    SECTION("Random density") {
        grid.setDensitiesRandom(1234);
        const LinearPiezoelectricSolver solver(grid, material);

        auto options = SolverOptions::defaults();
        options.maxIterations = 2000;

        auto results = solver.solve(options);
        const auto &cBar = results.cBar;
        const auto &epsilonBar = results.epsilonBar;
        const auto &dBar = results.dBar;

        REQUIRE(isSymmetric(cBar));
        REQUIRE(isPD(cBar));

        REQUIRE(isSymmetric(epsilonBar));
        REQUIRE(isPD(epsilonBar));

        SECTION("Voigt/Reuss bounds") {
            SECTION("cBar") {
                const auto voigt = elasticMaterial.voigt(grid);
                const auto reuss = elasticMaterial.reuss(grid);

                REQUIRE(reuss.trace() <= cBar.trace());
                REQUIRE(cBar.trace() <= voigt.trace());
            }

            SECTION("epsilonBar") {
                const auto voigt = dielectricMaterial.voigt(grid);
                const auto reuss = dielectricMaterial.reuss(grid);

                REQUIRE(reuss.trace() <= epsilonBar.trace());
                REQUIRE(epsilonBar.trace() <= voigt.trace());
            }
        }

        SECTION("Translational invariance") {
            Resolution shift;
            shift.fill(1);

            grid.translate(shift);

            const LinearPiezoelectricSolver solver2(grid, material);

            results = solver.solve();
            const auto &cBar2 = results.cBar;
            const auto &epsilonBar2 = results.epsilonBar;
            const auto &dBar2 = results.dBar;

            REQUIRE(cBar2.isApprox(cBar, NUMERICAL_ZERO));
            REQUIRE(epsilonBar2.isApprox(epsilonBar, NUMERICAL_ZERO));
            REQUIRE(dBar2.isApprox(dBar, NUMERICAL_ZERO));
        }
    }

    SECTION("Bad options") {
        grid.setDensitiesRandom(1234);
        const LinearPiezoelectricSolver solver(grid, material);

        SolverOptions options;
        options.maxIterations = 1;

        REQUIRE_THROWS_AS(solver.solve(options), std::runtime_error);
    }
}
