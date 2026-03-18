#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/material/material_aliases.hpp"
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

// Negative sign convention: linear dielectric material
// Positive sign convention: linear mass diffusive material
using Types = std::tuple<
    TypePair<Quad4Grid, LinearDielectricMaterial2d, DensityField2d, LinearDielectricSolver<Quad4Grid>>,
    TypePair<Quad4Grid, LinearMassDiffusiveMaterial2d, DensityField2d, LinearMassDiffusiveSolver<Quad4Grid>>,
    TypePair<Quad8Grid, LinearDielectricMaterial2d, DensityField2d, LinearDielectricSolver<Quad8Grid>>,
    TypePair<Quad8Grid, LinearMassDiffusiveMaterial2d, DensityField2d, LinearMassDiffusiveSolver<Quad8Grid>>,
    TypePair<Hex8Grid, LinearDielectricMaterial3d, DensityField3d, LinearDielectricSolver<Hex8Grid>>,
    TypePair<Hex8Grid, LinearMassDiffusiveMaterial3d, DensityField3d, LinearMassDiffusiveSolver<Hex8Grid>>,
    TypePair<Hex20Grid, LinearDielectricMaterial3d, DensityField3d, LinearDielectricSolver<Hex20Grid>>,
    TypePair<Hex20Grid, LinearMassDiffusiveMaterial3d, DensityField3d, LinearMassDiffusiveSolver<Hex20Grid>>
>;

TEMPLATE_LIST_TEST_CASE("monad::Solver: Test solve", "[monad]", Types) {
    using Grid = typename TestType::Grid;
    using Material = typename TestType::Material;
    using DensityField = typename TestType::DensityField;
    using Solver = typename TestType::Solver;

    using Resolution = typename Grid::Resolution;
    using Size = typename Grid::Size;
    using MaterialTensor = typename Material::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    Grid grid(resolution, size);

    MaterialTensor K = MaterialTensor::Random();
    // Make PSD
    K = K.transpose() * K;
    // Make PD
    K += MaterialTensor::Identity();

    const Material material(K);

    DensityField densityField(resolution);

    const Solver solver;

    SECTION("Density=1 → KBar=material") {
        densityField.setOnes();

        const auto results = solver.solve(grid, densityField, material);
        const auto &KBar = results.KBar;

        REQUIRE(KBar.isApprox(K, NUMERICAL_ZERO));
    }

    SECTION("Density=0 → KBar=0") {
        densityField.setZeros();

        const auto results = solver.solve(grid, densityField, material);
        const auto &KBar = results.KBar;

        // Set zero threshold to one order of magnitude larger
        REQUIRE(KBar.isZero(10.0 * NUMERICAL_ZERO));
    }

    SECTION("Random density") {
        densityField.setRandom(1234);

        const auto results = solver.solve(grid, densityField, material);
        const auto &KBar = results.KBar;

        REQUIRE(isSymmetric(KBar));
        REQUIRE(isPD(KBar));

        SECTION("Voigt/Reuss bounds") {
            const auto voigt = voigtBound(material, densityField);
            const auto reuss = reussBound(material, densityField);

            REQUIRE(reuss.trace() <= KBar.trace());
            REQUIRE(KBar.trace() <= voigt.trace());
        }

        SECTION("Translational invariance") {
            Resolution shift;
            shift.fill(1);

            densityField.translate(shift);

            const auto resultsShift = solver.solve(grid, densityField, material);

            const auto &KBarShift = resultsShift.KBar;

            REQUIRE(KBarShift.isApprox(KBar, NUMERICAL_ZERO));
        }
    }

    SECTION("Bad options") {
        densityField.setRandom(1234);

        SolverOptions options;
        options.maxIterations = 1;

        REQUIRE_THROWS_AS(solver.solve(grid, densityField, material, options), std::runtime_error);
    }
}
