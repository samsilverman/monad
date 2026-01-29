#include <tuple>
#include <stdexcept>
#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/fem/element/quad4.hpp"
#include "monad/fem/element/quad8.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::mechanical;
using namespace monad::detail;

using Types = std::tuple<Quad4, Quad8>;

TEMPLATE_LIST_TEST_CASE("monad::fem::mechanical::LinearElasticKernel (2d): Test bMatrix", "[monad]", Types) {
    using FieldMatrix = typename LinearElasticKernel<TestType>::FieldMatrix;

    auto nodes = TestType::localNodes();
    const auto points = TestType::quadratureRule().points;

    SECTION("BU=ε for unit strains") {
        const auto expected = Eigen::Matrix3d::Identity();

        for (auto point : points) {
            const auto B = LinearElasticKernel<TestType>::bMatrix(point, nodes);

            FieldMatrix U;
            U.setZero();

            // Displacements of unit axial strain ε₁₁
            U(Eigen::seq(0, Eigen::indexing::last, 2), 0) = nodes.col(0);

            // Displacements of unit axial strain ε₂₂
            U(Eigen::seq(1, Eigen::indexing::last, 2), 1) = nodes.col(1);

            // Displacements of unit shear strain ε₁₂
            U(Eigen::seq(0, Eigen::indexing::last, 2), 2) = 0.5 * nodes.col(1);
            U(Eigen::seq(1, Eigen::indexing::last, 2), 2) = 0.5 * nodes.col(0);

            REQUIRE((B * U).isApprox(expected, NUMERICAL_ZERO));
        }
    }

    SECTION("BU=0 for rigid body transformations") {
        for (auto point : points) {
            const auto B = LinearElasticKernel<TestType>::bMatrix(point, nodes);

            FieldMatrix U;
            U.setZero();

            // Unit x translation
            U(Eigen::seq(0, Eigen::indexing::last, 2), 0).setOnes();

            // Unit y translation
            U(Eigen::seq(1, Eigen::indexing::last, 2), 1).setOnes();

            // 90° ccw rotation
            U(Eigen::seq(0, Eigen::indexing::last, 2), 2) = -nodes.col(1);
            U(Eigen::seq(1, Eigen::indexing::last, 2), 2) = nodes.col(0);

            REQUIRE((B * U).isZero(NUMERICAL_ZERO));
        }
    }

    SECTION("Invalid element - degenerate element") {
        nodes = 0.0 * nodes;

        const auto point = nodes.row(0);

        REQUIRE_THROWS_AS(LinearElasticKernel<TestType>::bMatrix(point, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        const auto point = nodes.row(0);

        REQUIRE_THROWS_AS(LinearElasticKernel<TestType>::bMatrix(point, nodes), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::mechanical::LinearElasticKernel (2d): Test lhs/rhs", "[monad]", Types) {
    using FieldMatrix = typename LinearElasticKernel<TestType>::FieldMatrix;

    const LinearElasticMaterial2d material(1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    auto nodes = TestType::localNodes();

    SECTION("No errors") {
        const auto K = LinearElasticKernel<TestType>::lhs(material, nodes);

        REQUIRE(isSymmetric(K));
        REQUIRE(isPSD(K));

        SECTION("UᵀKU=UᵀF=0 for rigid body transformations") {
            const FieldMatrix F = LinearElasticKernel<TestType>::rhs(material, nodes);

            FieldMatrix U;
            U.setZero();

            // Unit x translation
            U(Eigen::seq(0, Eigen::indexing::last, 2), 0).setOnes();

            // Unit y translation
            U(Eigen::seq(1, Eigen::indexing::last, 2), 1).setOnes();

            // 90° ccw rotation
            U(Eigen::seq(0, Eigen::indexing::last, 2), 2) = -nodes.col(1);
            U(Eigen::seq(1, Eigen::indexing::last, 2), 2) = nodes.col(0);

            REQUIRE((U.transpose() * K * U).isZero(NUMERICAL_ZERO));
            REQUIRE((U.transpose() * F).isZero(NUMERICAL_ZERO));
        }
    }

    SECTION("Invalid element - degenerate element") {
        nodes = 0.0 * nodes;

        REQUIRE_THROWS_AS(LinearElasticKernel<TestType>::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearElasticKernel<TestType>::rhs(material, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        REQUIRE_THROWS_AS(LinearElasticKernel<TestType>::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearElasticKernel<TestType>::rhs(material, nodes), std::invalid_argument);
    }
}
