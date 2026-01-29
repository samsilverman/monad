#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include "monad/material/mechanical/linear_elastic_material_3d.hpp"
#include "monad/material/transport/linear_transport_material_aliases.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test initalization", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;

    SECTION("Invalid piezoelectric tensor") {
        const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
        const LinearDielectricMaterial3d dielectricMaterial(2.1);

        const CouplingTensor d {
            {0.0, 0.0, 0.0, 10.0, 0.0, 0.0},
            {0.0, 0.0, 10.0, 0.0, 0.0, 0.0},
            {10.0, 10.0, 10.0, 0.0, 0.0, 0.0}
        };

        REQUIRE_THROWS_AS(LinearPiezoelectricMaterial3d(elasticMaterial, dielectricMaterial, d), std::invalid_argument);
    }
}

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test elasticMaterial", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;

    const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
    const LinearDielectricMaterial3d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
        {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
        {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
    };

    const LinearPiezoelectricMaterial3d material(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material.elasticMaterial() == elasticMaterial);
}

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test dielectricMaterial", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;

    const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
    const LinearDielectricMaterial3d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
        {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
        {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
    };

    const LinearPiezoelectricMaterial3d material(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material.dielectricMaterial() == dielectricMaterial);
}

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test couplingTensor", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;

    const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
    const LinearDielectricMaterial3d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
        {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
        {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
    };

    const LinearPiezoelectricMaterial3d material(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material.couplingTensor().isApprox(d, NUMERICAL_ZERO));
}

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test materialTensor", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;
    using MaterialTensor = typename LinearPiezoelectricMaterial3d::MaterialTensor;

    const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
    const LinearDielectricMaterial3d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
        {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
        {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
    };

    const LinearPiezoelectricMaterial3d material(elasticMaterial, dielectricMaterial, d);

    MaterialTensor expected;
    expected << elasticMaterial.materialTensor(), -d.transpose(),
                -d, -dielectricMaterial.materialTensor();

    REQUIRE(material.materialTensor().isApprox(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test operator==", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;

    const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
    const LinearDielectricMaterial3d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
        {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
        {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
    };

    const LinearPiezoelectricMaterial3d material1(elasticMaterial, dielectricMaterial, d);
    const LinearPiezoelectricMaterial3d material2(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material1 == material2);
}

TEST_CASE("monad::LinearPiezoelectricMaterial3d: Test operator!=", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial3d::CouplingTensor;

    SECTION("Different elastic materials") {
        const LinearElasticMaterial3d elasticMaterial1(1.0, 0.3);
        const LinearElasticMaterial3d elasticMaterial2(1.1, 0.3);

        const LinearDielectricMaterial3d dielectricMaterial(2.1);

        const CouplingTensor d {
            {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
            {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
            {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
        };

        const LinearPiezoelectricMaterial3d material1(elasticMaterial1, dielectricMaterial, d);
        const LinearPiezoelectricMaterial3d material2(elasticMaterial2, dielectricMaterial, d);

        REQUIRE(material1 != material2);
    }

    SECTION("Different electric materials") {
        const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
        const LinearDielectricMaterial3d dielectricMaterial1(2.1);
        const LinearDielectricMaterial3d dielectricMaterial2(2.2);

        const CouplingTensor d {
            {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
            {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
            {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
        };

        const LinearPiezoelectricMaterial3d material1(elasticMaterial, dielectricMaterial1, d);
        const LinearPiezoelectricMaterial3d material2(elasticMaterial, dielectricMaterial2, d);

        REQUIRE(material1 != material2);
    }

    SECTION("Different piezoelectric tensors") {
        const LinearElasticMaterial3d elasticMaterial(1.0, 0.3);
        const LinearDielectricMaterial3d dielectricMaterial(2.1);

        const CouplingTensor d1 {
            {0.0, 0.0, 0.0, 0.01, 0.0, 0.0},
            {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
            {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
        };

        const CouplingTensor d2 {
            {0.0, 0.0, 0.0, 0.02, 0.0, 0.0},
            {0.0, 0.0, 0.01, 0.0, 0.0, 0.0},
            {0.01, 0.01, 0.01, 0.0, 0.0, 0.0}
        };

        const LinearPiezoelectricMaterial3d material1(elasticMaterial, dielectricMaterial, d1);
        const LinearPiezoelectricMaterial3d material2(elasticMaterial, dielectricMaterial, d2);

        REQUIRE(material1 != material2);
    }
}
