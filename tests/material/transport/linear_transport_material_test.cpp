#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/grid/quad4_grid.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

using Types = std::tuple<LinearTransportMaterial<2>, LinearTransportMaterial<3>>;

TEMPLATE_LIST_TEST_CASE("monad::LinearTransportMaterial: Test initalization - isotropic", "[monad]", Types) {
    SECTION("Invalid transport constant") {
        REQUIRE_THROWS_AS(TestType(0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType(-1.0), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::LinearTransportMaterial: Test initalization - anisotropic", "[monad]", Types) {
    using MaterialTensor = typename TestType::MaterialTensor;

    SECTION("Not Symmetric") {
        MaterialTensor K = MaterialTensor::Identity();
        K(0, 1) = 1;
 
        REQUIRE_THROWS_AS(TestType(K), std::invalid_argument);
    }

    SECTION("Symmetric but not PD") {
        MaterialTensor K = MaterialTensor::Identity();
        K(0, 0) = -5;

        REQUIRE_THROWS_AS(TestType(K), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::LinearTransportMaterial: Test materialTensor", "[monad]", Types) {
    using MaterialTensor = typename TestType::MaterialTensor;

    MaterialTensor K = MaterialTensor::Random();
    // Make PD
    K.noalias() = K * K.transpose() + MaterialTensor::Identity();
    
    const TestType material(K);

    REQUIRE(material.materialTensor().isApprox(K, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::LinearTransportMaterial: Test voigt/reuss", "[monad]", Types) {        
    const TestType material(2.1);

    // Just using to get density values (does not matter 2D vs. 3D)
    Quad4Grid grid({3, 3}, {1.0, 1.0});

    SECTION("Solid material: reuss=voigt=K") {
        const auto K = material.materialTensor();

        grid.setDensitiesOnes();

        const auto voigt = material.voigt(grid);
        const auto reuss = material.reuss(grid);

        REQUIRE(voigt.trace() == reuss.trace());
        REQUIRE(voigt.trace() == K.trace());
    }

    SECTION("Random material: reuss≤voigt") {
        grid.setDensitiesRandom(1234);
        
        const auto voigt = material.voigt(grid);
        const auto reuss = material.reuss(grid);        

        REQUIRE(reuss.trace() <= voigt.trace());
    }
}

TEMPLATE_LIST_TEST_CASE("monad::LinearTransportMaterial: Test operator==", "[monad]", Types) {
    const TestType material1(1.0);
    const TestType material2(1.0);

    REQUIRE(material1 == material2);
}

TEMPLATE_LIST_TEST_CASE("monad::LinearTransportMaterial: Test operator!=", "[monad]", Types) {
    const TestType material1(1.0);
    const TestType material2(1.1);

    REQUIRE(material1 != material2);
}
