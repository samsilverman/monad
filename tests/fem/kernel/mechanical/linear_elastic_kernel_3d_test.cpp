#include <tuple>
#include <stdexcept>
#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/fem/element/hex8.hpp"
#include "monad/fem/element/hex20.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/material/mechanical/linear_elastic_material_3d.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::mechanical;
using namespace monad::detail;

using Types = std::tuple<Hex8, Hex20>;

TEMPLATE_LIST_TEST_CASE("monad::fem::mechanical::LinearElasticKernel (3d): Test bMatrix", "[monad]", Types) {
    using FieldMatrix = typename LinearElasticKernel<TestType>::FieldMatrix;

    auto nodes = TestType::localNodes();
    const auto points = TestType::quadratureRule().points;

    SECTION("BU=ε for unit strains") {
        const auto expected = Eigen::Matrix<double, 6, 6>::Identity();

        for (auto point : points) {
            const auto B = LinearElasticKernel<TestType>::bMatrix(point, nodes);

            FieldMatrix U;
            U.setZero();

            // Displacements of unit axial strain ε₁₁
            U(Eigen::seq(0, Eigen::indexing::last, 3), 0) = nodes.col(0);

            // Displacements of unit axial strain ε₂₂
            U(Eigen::seq(1, Eigen::indexing::last, 3), 1) = nodes.col(1);

            // Displacements of unit axial strain ε₃₃
            U(Eigen::seq(2, Eigen::indexing::last, 3), 2) = nodes.col(2);

            // Displacements of unit shear strain ε₁₂
            U(Eigen::seq(0, Eigen::indexing::last, 3), 3) = 0.5 * nodes.col(1);
            U(Eigen::seq(1, Eigen::indexing::last, 3), 3) = 0.5 * nodes.col(0);

            // Displacements of unit shear strain ε₁₃
            U(Eigen::seq(0, Eigen::indexing::last, 3), 4) = 0.5 * nodes.col(2);
            U(Eigen::seq(2, Eigen::indexing::last, 3), 4) = 0.5 * nodes.col(0);

            // Displacements of unit shear strain ε₂₃
            U(Eigen::seq(1, Eigen::indexing::last, 3), 5) = 0.5 * nodes.col(2);
            U(Eigen::seq(2, Eigen::indexing::last, 3), 5) = 0.5 * nodes.col(1);

            REQUIRE((B * U).isApprox(expected, NUMERICAL_ZERO));
        }
    }

    SECTION("BU=0 for rigid body transformations") {
        for (auto point : points) {
            auto B = LinearElasticKernel<TestType>::bMatrix(point, nodes);

            FieldMatrix U;
            U.setZero();

            // Unit x translation
            U(Eigen::seq(0, Eigen::indexing::last, 3), 0).setOnes();

            // Unit y translation
            U(Eigen::seq(1, Eigen::indexing::last, 3), 1).setOnes();

            // Unit z translation
            U(Eigen::seq(2, Eigen::indexing::last, 3), 2).setOnes();

            // 90° ccw rotation in xy-plane
            U(Eigen::seq(0, Eigen::indexing::last, 3), 3) = -nodes.col(1);
            U(Eigen::seq(1, Eigen::indexing::last, 3), 3) = nodes.col(0);

            // 90° ccw rotation in xz-plane
            U(Eigen::seq(0, Eigen::indexing::last, 3), 4) = nodes.col(2);
            U(Eigen::seq(2, Eigen::indexing::last, 3), 4) = -nodes.col(0);

            // 90° ccw rotation in yz-plane
            U(Eigen::seq(1, Eigen::indexing::last, 3), 5) = -nodes.col(2);
            U(Eigen::seq(2, Eigen::indexing::last, 3), 5) = nodes.col(1);

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

TEMPLATE_LIST_TEST_CASE("monad::fem::mechanical::LinearElasticKernel (3d): Test lhs", "[monad]", Types) {
    using FieldMatrix = typename LinearElasticKernel<TestType>::FieldMatrix;

    const LinearElasticMaterial3d material(1.0, 0.3);
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
            U(Eigen::seq(0, Eigen::indexing::last, 3), 0).setOnes();

            // Unit y translation
            U(Eigen::seq(1, Eigen::indexing::last, 3), 1).setOnes();

            // Unit z translation
            U(Eigen::seq(2, Eigen::indexing::last, 3), 2).setOnes();

            // 90° ccw rotation in xy-plane
            U(Eigen::seq(0, Eigen::indexing::last, 3), 3) = -nodes.col(1);
            U(Eigen::seq(1, Eigen::indexing::last, 3), 3) = nodes.col(0);

            // 90° ccw rotation in xz-plane
            U(Eigen::seq(0, Eigen::indexing::last, 3), 4) = nodes.col(2);
            U(Eigen::seq(2, Eigen::indexing::last, 3), 4) = -nodes.col(0);

            // 90° ccw rotation in yz-plane
            U(Eigen::seq(1, Eigen::indexing::last, 3), 5) = -nodes.col(2);
            U(Eigen::seq(2, Eigen::indexing::last, 3), 5) = nodes.col(1);

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
