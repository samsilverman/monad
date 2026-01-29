#include <tuple>
#include <stdexcept>
#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/fem/element/quad4.hpp"
#include "monad/fem/element/quad8.hpp"
#include "monad/fem/kernel/multiphysics/linear_piezoelectric_kernel.hpp"
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/transport/linear_transport_material_aliases.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::mechanical;
using namespace monad::fem::scalar;
using namespace monad::fem::multiphysics;
using namespace monad::detail;

using Types = std::tuple<Quad4, Quad8>;

TEMPLATE_LIST_TEST_CASE("monad::fem::electrical::LinearPiezoelectricKernel (2d): Test lhs/rhs", "[monad]", Types) {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const LinearDielectricMaterial2d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.01, 0.0, 0.0},
        {0.0, 0.01, 0.01}
    };

    const LinearPiezoelectricMaterial2d material(elasticMaterial, dielectricMaterial, d);

    auto nodes = TestType::localNodes();

    SECTION("No errors") {
        const auto K = LinearPiezoelectricKernel<TestType>::lhs(material, nodes);

        REQUIRE(isSymmetric(K));
        REQUIRE(!isPSD(K));

        SECTION("xᵀKx=xᵀF=0 for rigid body transformations") {
            using FieldVector = Eigen::Vector<double, LinearPiezoelectricKernel<TestType>::NumDofs>;
            using MechanicalKernel = typename LinearPiezoelectricKernel<TestType>::MechanicalKernel;

            const auto F = LinearPiezoelectricKernel<TestType>::rhs(material, nodes);

            FieldVector x;

            const int numMechanicalDofs = MechanicalKernel::NumDofs;

            // Unit x translation
            x.setZero();
            x(Eigen::seq(0, numMechanicalDofs - 1, 2)).setOnes();

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // Unit y translation
            x.setZero();
            x(Eigen::seq(1, numMechanicalDofs - 1, 2)).setOnes();

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // 90° ccw rotation
            x.setZero();
            x(Eigen::seq(0, numMechanicalDofs - 1, 2)) = -nodes.col(1);
            x(Eigen::seq(1, numMechanicalDofs - 1, 2)) = nodes.col(0);

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // Unit electrical translation
            x.setZero();
            x(Eigen::seq(numMechanicalDofs, Eigen::indexing::last)).setOnes();

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));
        }
    }

    SECTION("Invalid element - degenerate element") {
        nodes = 0.0 * nodes;

        REQUIRE_THROWS_AS(LinearPiezoelectricKernel<TestType>::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearPiezoelectricKernel<TestType>::rhs(material, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        REQUIRE_THROWS_AS(LinearPiezoelectricKernel<TestType>::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearPiezoelectricKernel<TestType>::rhs(material, nodes), std::invalid_argument);
    }
}
