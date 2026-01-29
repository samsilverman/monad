#include <tuple>
#include <stdexcept>
#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/fem/element/quad4.hpp"
#include "monad/fem/element/quad8.hpp"
#include "monad/fem/element/hex8.hpp"
#include "monad/fem/element/hex20.hpp"
#include "monad/fem/kernel/scalar/linear_scalar_diffusive_kernel.hpp"
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::scalar;
using namespace monad::detail;

using Types = std::tuple<
    LinearScalarDiffusiveKernel<Quad4, GradientConvention::Negative>,
    LinearScalarDiffusiveKernel<Quad4, GradientConvention::Positive>,
    LinearScalarDiffusiveKernel<Quad8, GradientConvention::Negative>,
    LinearScalarDiffusiveKernel<Quad8, GradientConvention::Positive>,
    LinearScalarDiffusiveKernel<Hex8, GradientConvention::Negative>,
    LinearScalarDiffusiveKernel<Hex8, GradientConvention::Positive>,
    LinearScalarDiffusiveKernel<Hex20, GradientConvention::Negative>,
    LinearScalarDiffusiveKernel<Hex20, GradientConvention::Positive>
>;

TEMPLATE_LIST_TEST_CASE("monad::fem::scalar::LinearScalarDiffusiveKernel: Test bMatrix", "[monad]", Types) {
    using Element = typename TestType::Element;

    auto nodes = Element::localNodes();
    const auto points = Element::quadratureRule().points;

    SECTION("Bφ=∇φ for unit scalar potentials") {
        const auto expected = Eigen::Matrix<double, Element::Dim, Element::Dim>::Identity();

        for (auto point : points) {
            const auto B = TestType::bMatrix(point, nodes);

            // Scalar potential of unit axial scalar potential gradient field ∇φᵢᵢ
            const auto Phi = TestType::GradSign * nodes;

            REQUIRE((B * Phi).isApprox(expected, NUMERICAL_ZERO));
        }
    }

    SECTION("Bφ=0 for rigid body transformations") {
        using Element = typename TestType::Element;
        using FieldVector = Eigen::Vector<double, Element::NumNodes>;

        for (auto point : points) {
            const auto B = TestType::bMatrix(point, nodes);

            // Unit transformation
            const FieldVector phi = FieldVector::Ones();

            REQUIRE((B * phi).isZero(NUMERICAL_ZERO));
        }
    }

    SECTION("Invalid element - degenerate element") {
        nodes = 0.0 * nodes;

        const auto point = nodes.row(0);

        REQUIRE_THROWS_AS(TestType::bMatrix(point, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        const auto point = nodes.row(0);

        REQUIRE_THROWS_AS(TestType::bMatrix(point, nodes), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::electrical::LinearScalarDiffusiveKernel: Test lhs/rhs", "[monad]", Types) {
    using Element = typename TestType::Element;

    const LinearTransportMaterial<Element::Dim> material(2.1);
    auto nodes = Element::localNodes();

    SECTION("No errors") {
        const auto K = TestType::lhs(material, nodes);

        REQUIRE(isSymmetric(K));
        REQUIRE(isPSD(K));

        SECTION("φᵀKφ=φᵀF=0 for rigid body transformations") {
            using FieldVector = Eigen::Vector<double, Element::NumNodes>;

            const auto F = TestType::rhs(material, nodes);

            // Unit translation
            const FieldVector phi = FieldVector::Ones();

            REQUIRE((phi.transpose() * K * phi).isZero(NUMERICAL_ZERO));
            REQUIRE((phi.transpose() * F).isZero(NUMERICAL_ZERO));
        }
    }

    SECTION("Invalid element - degenerate element") {
        nodes = 0.0 * nodes;

        REQUIRE_THROWS_AS(TestType::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType::rhs(material, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        REQUIRE_THROWS_AS(TestType::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType::rhs(material, nodes), std::invalid_argument);
    }
}
