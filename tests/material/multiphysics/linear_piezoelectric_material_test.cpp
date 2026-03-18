#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/material_aliases.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

using Types = std::tuple<LinearPiezoelectricMaterial2d, LinearPiezoelectricMaterial3d>;

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test initalization", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;

    SECTION("Invalid stiffness tensor") {
        StiffnessTensor c = StiffnessTensor::Random();
        // Make PSD
        c = c.transpose() * c;
        c(0, 0) = 0.0;

        PermittivityTensor epsilon = PermittivityTensor::Random();
        // Make PSD
        epsilon = epsilon.transpose() * epsilon;
        // Make PD
        epsilon += PermittivityTensor::Identity();

        const CouplingTensor d = CouplingTensor::Constant(0.01);

        REQUIRE_THROWS_AS(TestType(c, epsilon, d), std::invalid_argument);
    }

    SECTION("Invalid permittivity tensor") {
        StiffnessTensor c = StiffnessTensor::Random();
        // Make PSD
        c = c.transpose() * c;
        // Make PD
        c += StiffnessTensor::Identity();

        PermittivityTensor epsilon = PermittivityTensor::Random();
        // Make PSD
        epsilon = epsilon.transpose() * epsilon;
        epsilon(0, 0) = 0.0;

        const CouplingTensor d = CouplingTensor::Constant(0.01);

        REQUIRE_THROWS_AS(TestType(c, epsilon, d), std::invalid_argument);
    }

    SECTION("Invalid piezoelectric coupling tensor") {
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

        const CouplingTensor d = CouplingTensor::Constant(10.0);

        REQUIRE_THROWS_AS(TestType(c, epsilon, d), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test stiffnessTensor", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;

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

    const TestType material(c, epsilon, d);

    REQUIRE(material.stiffnessTensor().isApprox(c, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test permittivityTensor", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;

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

    const TestType material(c, epsilon, d);

    REQUIRE(material.permittivityTensor().isApprox(epsilon, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test couplingTensor", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;

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

    const TestType material(c, epsilon, d);

    REQUIRE(material.couplingTensor().isApprox(d, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test materialTensor", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;
    using MaterialTensor = typename TestType::MaterialTensor;

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

    const TestType material(c, epsilon, d);

    MaterialTensor op;
    op << c, -d.transpose(),
          -d, -epsilon;

    REQUIRE(material.materialTensor().isApprox(op, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test operator==", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;

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

    const TestType material1(c, epsilon, d);
    const TestType material2(c, epsilon, d);

    REQUIRE(material1 == material2);
}

TEMPLATE_LIST_TEST_CASE("monad::LinearPiezoelectricMaterial: Test operator!=", "[monad]", Types) {
    using StiffnessTensor = typename TestType::StiffnessTensor;
    using PermittivityTensor = typename TestType::PermittivityTensor;
    using CouplingTensor = typename TestType::CouplingTensor;

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

    CouplingTensor d = CouplingTensor::Constant(0.01);

    const TestType material1(c, epsilon, d);

    SECTION("Different stiffness tensor") {
        c(1, 1) += 0.01;

        const TestType material2(c, epsilon, d);

        REQUIRE(material1 != material2);
    }
    
    SECTION("Different permittivity tensor") {
        epsilon(1, 1) += 0.01;

        const TestType material2(c, epsilon, d);

        REQUIRE(material1 != material2);
    }

    SECTION("Different piezoelectric coupling tensor") {
        d(1, 1) += 0.01;

        const TestType material2(c, epsilon, d);

        REQUIRE(material1 != material2);
    }
}
