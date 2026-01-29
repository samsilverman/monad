#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/transport/linear_transport_material_aliases.hpp"
#include "monad/material/multiphysics/linear_piezoelectric_material.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test initalization", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    SECTION("Invalid piezoelectric tensor") {
        const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
        const LinearDielectricMaterial2d dielectricMaterial(2.1);

        const CouplingTensor d {
            {10.0, 0.0, 0.0},
            {0.0, 10.0, 10.0}
        };

        REQUIRE_THROWS_AS(LinearPiezoelectricMaterial2d(elasticMaterial, dielectricMaterial, d), std::invalid_argument);
    }
}

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test elasticMaterial", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const LinearDielectricMaterial2d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.01, 0.0, 0.0},
        {0.0, 0.01, 0.01}
    };

    const LinearPiezoelectricMaterial2d material(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material.elasticMaterial() == elasticMaterial);
}

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test dielectricMaterial", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const LinearDielectricMaterial2d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.01, 0.0, 0.0},
        {0.0, 0.01, 0.01}
    };

    const LinearPiezoelectricMaterial2d material(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material.dielectricMaterial() == dielectricMaterial);
}

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test couplingTensor", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const LinearDielectricMaterial2d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.01, 0.0, 0.0},
        {0.0, 0.01, 0.01}
    };

    const LinearPiezoelectricMaterial2d material(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material.couplingTensor().isApprox(d, NUMERICAL_ZERO));
}

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test materialTensor", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;
    using MaterialTensor = typename LinearPiezoelectricMaterial2d::MaterialTensor;

    const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const LinearDielectricMaterial2d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.01, 0.0, 0.0},
        {0.0, 0.01, 0.01}
    };

    const LinearPiezoelectricMaterial2d material(elasticMaterial, dielectricMaterial, d);

    MaterialTensor expected;
    expected << elasticMaterial.materialTensor(), -d.transpose(),
                -d, -dielectricMaterial.materialTensor();

    REQUIRE(material.materialTensor().isApprox(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test operator==", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    const LinearDielectricMaterial2d dielectricMaterial(2.1);

    const CouplingTensor d {
        {0.01, 0.0, 0.0},
        {0.0, 0.01, 0.01}
    };

    const LinearPiezoelectricMaterial2d material1(elasticMaterial, dielectricMaterial, d);
    const LinearPiezoelectricMaterial2d material2(elasticMaterial, dielectricMaterial, d);

    REQUIRE(material1 == material2);
}

TEST_CASE("monad::LinearPiezoelectricMaterial2d: Test operator!=", "[monad]") {
    using CouplingTensor = typename LinearPiezoelectricMaterial2d::CouplingTensor;

    SECTION("Different elastic materials") {
        const LinearElasticMaterial2d elasticMaterial1(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
        const LinearElasticMaterial2d elasticMaterial2(1.1, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
        const LinearDielectricMaterial2d dielectricMaterial(2.1);

        const CouplingTensor d {
            {0.01, 0.0, 0.0},
            {0.0, 0.01, 0.01}
        };

        LinearPiezoelectricMaterial2d material1(elasticMaterial1, dielectricMaterial, d);
        LinearPiezoelectricMaterial2d material2(elasticMaterial2, dielectricMaterial, d);

        REQUIRE(material1 != material2);
    }

    SECTION("Different electric materials") {
        const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
        const LinearDielectricMaterial2d dielectricMaterial1(2.1);
        const LinearDielectricMaterial2d dielectricMaterial2(2.2);

        const CouplingTensor d {
            {0.01, 0.0, 0.0},
            {0.0, 0.01, 0.01}
        };

        const LinearPiezoelectricMaterial2d material1(elasticMaterial, dielectricMaterial1, d);
        const LinearPiezoelectricMaterial2d material2(elasticMaterial, dielectricMaterial2, d);

        REQUIRE(material1 != material2);
    }

    SECTION("Different piezoelectric tensors") {
        const LinearElasticMaterial2d elasticMaterial(1.0, 0.3, monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress);
        const LinearDielectricMaterial2d dielectricMaterial(2.1);

        const CouplingTensor d1 {
            {0.01, 0.0, 0.0},
            {0.0, 0.01, 0.01}
        };

        const CouplingTensor d2 {
            {0.02, 0.0, 0.0},
            {0.0, 0.01, 0.01}
        };

        const LinearPiezoelectricMaterial2d material1(elasticMaterial, dielectricMaterial, d1);
        const LinearPiezoelectricMaterial2d material2(elasticMaterial, dielectricMaterial, d2);

        REQUIRE(material1 != material2);
    }
}
