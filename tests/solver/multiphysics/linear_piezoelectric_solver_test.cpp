#include <tuple>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/mechanical/linear_elastic_material_3d.hpp"
#include "monad/material/material_aliases.hpp"
#include "monad/material/bounds.hpp"
#include "monad/field/field_aliases.hpp"
#include "monad/solver/solver_aliases.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::detail;

template <class GridT, class MechanicalMaterialT, class ElectricalMaterialT, class MaterialT, class DensityFieldT, class SolverT>
struct TypePair {
    using Grid = GridT;
    using MechanicalMaterial = MechanicalMaterialT;
    using ElectricalMaterial = ElectricalMaterialT;
    using Material = MaterialT;
    using DensityField = DensityFieldT;
    using Solver = SolverT;
};

using Types = std::tuple<
    TypePair<Quad4Grid, LinearElasticMaterial2d, LinearDielectricMaterial2d, LinearPiezoelectricMaterial2d, DensityField2d, LinearPiezoelectricSolver<Quad4Grid>>,
    TypePair<Quad8Grid, LinearElasticMaterial2d, LinearDielectricMaterial2d, LinearPiezoelectricMaterial2d, DensityField2d, LinearPiezoelectricSolver<Quad8Grid>>,
    TypePair<Hex8Grid, LinearElasticMaterial3d, LinearDielectricMaterial3d, LinearPiezoelectricMaterial3d, DensityField3d, LinearPiezoelectricSolver<Hex8Grid>>,
    TypePair<Hex20Grid, LinearElasticMaterial3d, LinearDielectricMaterial3d, LinearPiezoelectricMaterial3d, DensityField3d, LinearPiezoelectricSolver<Hex20Grid>>
>;

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricSolver: Test solve", "[monad]", Types) {
    using Grid = typename TestType::Grid;
    using MechanicalMaterial = typename TestType::MechanicalMaterial;
    using ElectricalMaterial = typename TestType::ElectricalMaterial;
    using Material = typename TestType::Material;
    using DensityField = typename TestType::DensityField;
    using Solver = typename TestType::Solver;

    using Resolution = typename Grid::Resolution;
    using Size = typename Grid::Size;
    using CouplingTensor = typename Material::CouplingTensor;
    using StiffnessTensor = typename MechanicalMaterial::MaterialTensor;
    using PermittivityTensor = typename ElectricalMaterial::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    Grid grid(resolution, size);

    StiffnessTensor c = StiffnessTensor::Random();
    // Make PSD
    c = c.transpose() * c;
    // Make PD
    c += StiffnessTensor::Identity();

    const MechanicalMaterial elasticMaterial(c);

    PermittivityTensor epsilon = PermittivityTensor::Random();
    // Make PSD
    epsilon = epsilon.transpose() * epsilon;
    // Make PD
    epsilon += PermittivityTensor::Identity();

    const ElectricalMaterial dielectricMaterial(epsilon);

    const CouplingTensor d = 0.1 * CouplingTensor::Random();

    const Material material(elasticMaterial, dielectricMaterial, d);

    DensityField densityField(resolution);

    const Solver solver;

    SECTION("Density=1 → homogenized=material") {
        densityField.setOnes();

        const auto results = solver.solve(grid, densityField, material);
        const auto &cBar = results.cBar;
        const auto &epsilonBar = results.epsilonBar;
        const auto &dBar = results.dBar;

        REQUIRE(cBar.isApprox(c, NUMERICAL_ZERO));
        REQUIRE(epsilonBar.isApprox(epsilon, NUMERICAL_ZERO));
        REQUIRE(dBar.isApprox(d, NUMERICAL_ZERO));
    }

    SECTION("Density=0 → CBar=0") {
        densityField.setZeros();

        const auto results = solver.solve(grid, densityField, material);
        const auto &cBar = results.cBar;
        const auto &epsilonBar = results.epsilonBar;
        const auto &dBar = results.dBar;

        // Set zero threshold to one order of magnitude larger
        REQUIRE(cBar.isZero(10.0 * NUMERICAL_ZERO));
        REQUIRE(epsilonBar.isZero(10.0 * NUMERICAL_ZERO));
        REQUIRE(dBar.isZero(10.0 * NUMERICAL_ZERO));
    }

    SECTION("Random density") {
        densityField.setRandom(1234);

        const auto results = solver.solve(grid, densityField, material);
        const auto &opBar = results.opBar;

        REQUIRE(isSymmetric(opBar));
        REQUIRE(!isPD(opBar));

        SECTION("Voigt/Reuss bounds") {
            const auto voigt = voigtBound(material, densityField);
            const auto reuss = reussBound(material, densityField);

            REQUIRE(reuss.trace() <= opBar.trace());
            REQUIRE(opBar.trace() <= voigt.trace());
        }

        SECTION("Translational invariance") {
            Resolution shift;
            shift.fill(1);

            densityField.translate(shift);

            const auto resultsShift = solver.solve(grid, densityField, material);

            const auto &opBarShift = resultsShift.opBar;

            REQUIRE(opBarShift.isApprox(opBar, NUMERICAL_ZERO));
        }
    }

    SECTION("Bad options") {
        densityField.setRandom(1234);

        SolverOptions options;
        options.maxIterations = 1;

        REQUIRE_THROWS_AS(solver.solve(grid, densityField, material, options), std::runtime_error);
    }
}
