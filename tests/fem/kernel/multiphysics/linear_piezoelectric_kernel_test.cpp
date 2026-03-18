#include <tuple>
#include <stdexcept>
#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/fem/element/hex8.hpp"
#include "monad/fem/element/hex20.hpp"
#include "monad/fem/element/quad4.hpp"
#include "monad/fem/element/quad8.hpp"
#include "monad/fem/kernel/multiphysics/linear_piezoelectric_kernel.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::mechanical;
using namespace monad::fem::scalar;
using namespace monad::fem::multiphysics;
using namespace monad::detail;
using namespace monad::material;

using Types2d = std::tuple<Quad4, Quad8>;
using Types3d = std::tuple<Hex8, Hex20>;

TEMPLATE_LIST_TEST_CASE("monad::fem::electrical::LinearPiezoelectricKernel (2d): Test lhs/rhs", "[monad]", Types2d) {
    using Kernel = LinearPiezoelectricKernel<TestType>;
    using Material = typename Kernel::Material;
    using StiffnessTensor = typename Material::StiffnessTensor;
    using PermittivityTensor = typename Material::PermittivityTensor;
    using CouplingTensor = typename Material::CouplingTensor;

    StiffnessTensor c = StiffnessTensor::Random();
    // Make PSD
    c = c.transpose() * c;
    // Make PD
    c += StiffnessTensor::Identity();

    PermittivityTensor epsilon = PermittivityTensor::Random();
    // Make PSD
    epsilon = epsilon.transpose() * epsilon;
    // Make PD
    epsilon += PermittivityTensor::Identity();

    const CouplingTensor d = CouplingTensor::Constant(0.01);

    const Material material(c, epsilon, d);

    auto nodes = TestType::localNodes();

    SECTION("No errors") {
        const auto K = Kernel::lhs(material, nodes);

        REQUIRE(isSymmetric(K));
        REQUIRE(!isPSD(K));

        SECTION("xᵀKx=xᵀF=0 for rigid body transformations") {
            using FieldVector = Eigen::Vector<double, Kernel::NumDofs>;
            using MechanicalKernel = typename Kernel::MechanicalKernel;

            const auto F = Kernel::rhs(material, nodes);

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

        REQUIRE_THROWS_AS(Kernel::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(Kernel::rhs(material, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        REQUIRE_THROWS_AS(Kernel::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(Kernel::rhs(material, nodes), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::fem::electrical::LinearPiezoelectricKernel (3d): Test lhs/rhs", "[monad]", Types3d) {
    using Kernel = LinearPiezoelectricKernel<TestType>;
    using Material = typename Kernel::Material;
    using StiffnessTensor = typename Material::StiffnessTensor;
    using PermittivityTensor = typename Material::PermittivityTensor;
    using CouplingTensor = typename Material::CouplingTensor;

    StiffnessTensor c = StiffnessTensor::Random();
    // Make PSD
    c = c.transpose() * c;
    // Make PD
    c += StiffnessTensor::Identity();

    PermittivityTensor epsilon = PermittivityTensor::Random();
    // Make PSD
    epsilon = epsilon.transpose() * epsilon;
    // Make PD
    epsilon += PermittivityTensor::Identity();

    const CouplingTensor d = CouplingTensor::Constant(0.01);

    const Material material(c, epsilon, d);

    auto nodes = TestType::localNodes();

    SECTION("No errors") {
        const auto K = Kernel::lhs(material, nodes);

        REQUIRE(isSymmetric(K));
        REQUIRE(!isPSD(K));

        SECTION("xᵀKx=xᵀF=0 for rigid body transformations") {
            using FieldVector = Eigen::Vector<double, Kernel::NumDofs>;
            using MechanicalKernel = typename Kernel::MechanicalKernel;

            const auto F = Kernel::rhs(material, nodes);

            FieldVector x;

            const int numMechanicalDofs = MechanicalKernel::NumDofs;

            // Unit x translation
            x.setZero();
            x(Eigen::seq(0, numMechanicalDofs - 1, 3)).setOnes();

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // Unit y translation
            x.setZero();
            x(Eigen::seq(1, numMechanicalDofs - 1, 3)).setOnes();

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // Unit z translation
            x.setZero();
            x(Eigen::seq(2, numMechanicalDofs - 1, 3)).setOnes();

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // 90° ccw rotation in xy-plane
            x.setZero();
            x(Eigen::seq(0, numMechanicalDofs - 1, 3)) = -nodes.col(1);
            x(Eigen::seq(1, numMechanicalDofs - 1, 3)) = nodes.col(0);

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // 90° ccw rotation in xz-plane
            x.setZero();
            x(Eigen::seq(0, numMechanicalDofs - 1, 3)) = nodes.col(2);
            x(Eigen::seq(2, numMechanicalDofs - 1, 3)) = -nodes.col(0);

            REQUIRE((x.transpose() * K * x).isZero(NUMERICAL_ZERO));
            REQUIRE((x.transpose() * F).isZero(NUMERICAL_ZERO));

            // 90° ccw rotation in yz-plane
            x.setZero();
            x(Eigen::seq(1, numMechanicalDofs - 1, 3)) = -nodes.col(2);
            x(Eigen::seq(2, numMechanicalDofs - 1, 3)) = nodes.col(1);

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

        REQUIRE_THROWS_AS(Kernel::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(Kernel::rhs(material, nodes), std::invalid_argument);
    }

    SECTION("Invalid element - reverse node ordering") {
        nodes = nodes.rowwise().reverse();

        REQUIRE_THROWS_AS(Kernel::lhs(material, nodes), std::invalid_argument);
        REQUIRE_THROWS_AS(Kernel::rhs(material, nodes), std::invalid_argument);
    }
}
