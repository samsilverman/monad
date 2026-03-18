#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/mechanical/linear_elastic_material_3d.hpp"
#include "monad/material/bounds.hpp"
#include "monad/field/field_aliases.hpp"
#include "monad/solver/solver_aliases.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::detail;

template <class GridT, class MaterialT, class DensityFieldT, class SolverT>
struct TypePair {
    using Grid = GridT;
    using Material = MaterialT;
    using DensityField = DensityFieldT;
    using Solver = SolverT;
};

using Types = std::tuple<
    TypePair<Quad4Grid, LinearElasticMaterial2d, DensityField2d, LinearElasticSolver<Quad4Grid>>,
    TypePair<Quad8Grid, LinearElasticMaterial2d, DensityField2d, LinearElasticSolver<Quad8Grid>>,
    TypePair<Hex8Grid, LinearElasticMaterial3d, DensityField3d, LinearElasticSolver<Hex8Grid>>,
    TypePair<Hex20Grid, LinearElasticMaterial3d, DensityField3d, LinearElasticSolver<Hex20Grid>>
>;

TEMPLATE_LIST_TEST_CASE("monad::LinearElasticSolver: Test solve", "[monad]", Types) {
    using Grid = typename TestType::Grid;
    using Material = typename TestType::Material;
    using DensityField = typename TestType::DensityField;
    using Solver = typename TestType::Solver;

    using Resolution = typename Grid::Resolution;
    using Size = typename Grid::Size;
    using StiffnessTensor = typename Material::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    Grid grid(resolution, size);
    
    StiffnessTensor C = StiffnessTensor::Random();
    // Make PSD
    C = C.transpose() * C;
    // Make PD
    C += StiffnessTensor::Identity();

    const Material material(C);

    DensityField densityField(resolution);

    const Solver solver;

    SECTION("Density=1 → CBar=C") {
        densityField.setOnes();

        const auto results = solver.solve(grid, densityField, material);
        const auto &CBar = results.CBar;

        REQUIRE(CBar.isApprox(C, NUMERICAL_ZERO));
    }

    SECTION("Density=0 → CBar=0") {
        densityField.setZeros();

        const auto results = solver.solve(grid, densityField, material);
        const auto &CBar = results.CBar;

        // Set zero threshold to one order of magnitude larger
        REQUIRE(CBar.isZero(10.0 * NUMERICAL_ZERO));
    }

    SECTION("Random density") {
        densityField.setRandom(1234);

        const auto results = solver.solve(grid, densityField, material);
        const auto &CBar = results.CBar;

        REQUIRE(isSymmetric(CBar));
        REQUIRE(isPD(CBar));

        SECTION("Voigt/Reuss bounds") {
            const auto voigt = voigtBound(material, densityField);
            const auto reuss = reussBound(material, densityField);

            REQUIRE(reuss.trace() <= CBar.trace());
            REQUIRE(CBar.trace() <= voigt.trace());
        }

        SECTION("Translational invariance") {
            Resolution shift;
            shift.fill(1);

            densityField.translate(shift);

            const auto resultsShift = solver.solve(grid, densityField, material);

            const auto &CBarShift = resultsShift.CBar;

            REQUIRE(CBarShift.isApprox(CBar, NUMERICAL_ZERO));
        }
    }

    SECTION("Bad options") {
        densityField.setRandom(1234);

        SolverOptions options;
        options.maxIterations = 1;

        REQUIRE_THROWS_AS(solver.solve(grid, densityField, material, options), std::runtime_error);
    }
}
