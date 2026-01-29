#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/grid/quad4_grid.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

TEST_CASE("monad::LinearElasticMaterial2d: Test initalization - isotropic", "[monad]") {
    SECTION("Invalid Young's modulus") {
        REQUIRE_THROWS_AS(LinearElasticMaterial2d(0.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearElasticMaterial2d(-1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress), std::invalid_argument);

        REQUIRE_THROWS_AS(LinearElasticMaterial2d(0.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStrain), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearElasticMaterial2d(-1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStrain), std::invalid_argument);
    }

    SECTION("Invalid Poisson's ratio") {
        REQUIRE_THROWS_AS(LinearElasticMaterial2d(1.0, 0.5, LinearElasticMaterial2d::PlaneCondition::PlaneStress), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearElasticMaterial2d(1.0, -1.0, LinearElasticMaterial2d::PlaneCondition::PlaneStress), std::invalid_argument);

        REQUIRE_THROWS_AS(LinearElasticMaterial2d(1.0, 0.5, LinearElasticMaterial2d::PlaneCondition::PlaneStrain), std::invalid_argument);
        REQUIRE_THROWS_AS(LinearElasticMaterial2d(1.0, -1.0, LinearElasticMaterial2d::PlaneCondition::PlaneStrain), std::invalid_argument);
    }
}

TEST_CASE("monad::LinearElasticMaterial2d: Test initalization - anisotropic", "[monad]") {
    using MaterialTensor = typename LinearElasticMaterial2d::MaterialTensor;

    SECTION("Not Symmetric") {
        MaterialTensor C = MaterialTensor::Identity();
        C(0, 1) = 1;

        REQUIRE_THROWS_AS(LinearElasticMaterial2d(C), std::invalid_argument);
    }

    SECTION("Symmetric but not PD") {
        MaterialTensor C = MaterialTensor::Identity();
        C(0, 0) = -5;

        REQUIRE_THROWS_AS(LinearElasticMaterial2d(C), std::invalid_argument);
    }
}

TEST_CASE("monad::LinearElasticMaterial2d: Test materialTensor", "[monad]") {
    using MaterialTensor = typename LinearElasticMaterial2d::MaterialTensor;

    MaterialTensor C = MaterialTensor::Random();
    // Make PD
    C.noalias() = C * C.transpose() + MaterialTensor::Identity();

    const LinearElasticMaterial2d material(C);

    REQUIRE(material.materialTensor().isApprox(C, NUMERICAL_ZERO));
}

TEST_CASE("monad::LinearElasticMaterial2d: Test voigt/reuss", "[monad]") {        
    const LinearElasticMaterial2d material(1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress);

    // Just using to get density values (does not matter 2D vs. 3D)
    Quad4Grid grid({3, 3}, {1.0, 1.0});

    SECTION("Solid material: reuss=voigt=C") {
        const auto C = material.materialTensor();

        grid.setDensitiesOnes();
        
        const auto voigt = material.voigt(grid);
        const auto reuss = material.reuss(grid);
        
        REQUIRE(voigt.trace() == reuss.trace());
        REQUIRE(voigt.trace() == C.trace());
    }

    SECTION("Random material: reuss≤voigt") {
        grid.setDensitiesRandom(1234);
        
        const auto voigt = material.voigt(grid);
        const auto reuss = material.reuss(grid);        

        REQUIRE(reuss.trace() <= voigt.trace());
    }
}

TEST_CASE("monad::LinearElasticMaterial2d: Test operator==", "[monad]") {
    LinearElasticMaterial2d material1(1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    LinearElasticMaterial2d material2(1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress);

    REQUIRE(material1 == material2);
}

TEST_CASE("monad::LinearElasticMaterial2d: Test operator!=", "[monad]") {
    LinearElasticMaterial2d material1(1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    LinearElasticMaterial2d material2(1.1, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    LinearElasticMaterial2d material3(1.0, 0.4, LinearElasticMaterial2d::PlaneCondition::PlaneStress);
    LinearElasticMaterial2d material4(1.0, 0.3, LinearElasticMaterial2d::PlaneCondition::PlaneStrain);

    REQUIRE(material1 != material2);
    REQUIRE(material1 != material3);
    REQUIRE(material1 != material4);
}
