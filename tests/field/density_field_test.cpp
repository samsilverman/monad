#include <tuple>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "monad/field/field_aliases.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

using Types = std::tuple<DensityField2d, DensityField3d>;

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test initalization", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    SECTION("Invalid resolution") {
        for (std::size_t i = 0; i < TestType::Dim; ++i) {
            Resolution resolutionInvalid = resolution;
            resolutionInvalid[i] = 0;

            REQUIRE_THROWS_AS(TestType(resolutionInvalid), std::invalid_argument);
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test resolution", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    const TestType densityField(resolution);

    const Resolution actual = densityField.resolution();

    REQUIRE(actual == resolution);
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test numElements", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    const TestType densityField(resolution);

    const std::size_t actual = static_cast<std::size_t>(std::pow(2.0, TestType::Dim));
    const std::size_t expected = densityField.numElements();

    REQUIRE(actual == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test densities", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using DensityList = typename TestType::DensityList;

    Resolution resolution;
    resolution.fill(2);

    const TestType densityField(resolution);

    const DensityList actual(densityField.numElements(), NUMERICAL_ZERO);
    const DensityList expected = densityField.densities();

    REQUIRE_THAT(actual, Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test getDensity/setDensity", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField(resolution);

    SECTION("No errors") {
        REQUIRE_THAT(densityField.getDensity(1), Catch::Matchers::WithinAbs(NUMERICAL_ZERO, NUMERICAL_ZERO));

        densityField.setDensity(1, 0.4);

        REQUIRE_THAT(densityField.getDensity(1), Catch::Matchers::WithinAbs(0.4, NUMERICAL_ZERO));
    }

    SECTION("Invalid index") {
        const std::size_t badIndex = densityField.numElements();

        REQUIRE_THROWS_AS(densityField.getDensity(badIndex), std::out_of_range);
        REQUIRE_THROWS_AS(densityField.setDensity(badIndex, 0.4), std::out_of_range);
    }

    SECTION("Invalid density") {
        REQUIRE_THROWS_AS(densityField.setDensity(0, -0.1), std::invalid_argument);
        REQUIRE_THROWS_AS(densityField.setDensity(0, 1.1), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test setDensities", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using DensityList = typename TestType::DensityList;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField(resolution);

    DensityList densities(densityField.numElements(), 0.4);

    SECTION("No errors") {
        densityField.setDensities(densities);

        REQUIRE_THAT(densityField.densities(), Catch::Matchers::Approx(densities).margin(NUMERICAL_ZERO));
    }

    SECTION("Invalid density") {
        densities[0] = -0.1;
        
        REQUIRE_THROWS_AS(densityField.setDensities(densities), std::invalid_argument);
        
        densities[0] = 1.1;
        
        REQUIRE_THROWS_AS(densityField.setDensities(densities), std::invalid_argument);
    }

    SECTION("Invalid size") {
        densities.push_back(0.4);

        REQUIRE_THROWS_AS(densityField.setDensities(densities), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test setConstant", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using DensityList = typename TestType::DensityList;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField(resolution);

    SECTION("No errors") {
        DensityList expected(densityField.numElements(), 0.4);

        densityField.setConstant(0.4);

        REQUIRE_THAT(densityField.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
    }

    SECTION("Invalid density") {        
        REQUIRE_THROWS_AS(densityField.setConstant(-0.1), std::invalid_argument);  
        REQUIRE_THROWS_AS(densityField.setConstant(1.1), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test setZeros", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField1(resolution);
    TestType densityField2(resolution);

    densityField1.setConstant(0.0);
    densityField2.setZeros();

    REQUIRE(densityField1 == densityField2);
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test setOnes", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField1(resolution);
    TestType densityField2(resolution);

    densityField1.setConstant(1.0);
    densityField2.setOnes();

    REQUIRE(densityField1 == densityField2);
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test setRandom", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using DensityList = typename TestType::DensityList;
    
    Resolution resolution;
    resolution.fill(2);

    TestType densityField(resolution);

    densityField.setRandom(1234);

    if constexpr (TestType::Dim == 2) {
        DensityList expected{
            0.19151945020806033,
            0.49766366637723142,
            0.62210876648829061,
            0.8178384429351051,
        };

        REQUIRE_THAT(densityField.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
    }

    else {
        DensityList expected{
            0.19151945020806033,
            0.49766366637723142,
            0.62210876648829061,
            0.8178384429351051,
            0.43772773734240972,
            0.61211189362502472,
            0.78535858373748102,
            0.77135991905149071
        };

        REQUIRE_THAT(densityField.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
    }
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test translate", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    
    Resolution resolution;
    resolution.fill(2);

    TestType densityField(resolution);

    densityField.setDensity(1, 0.4);

    if constexpr (TestType::Dim == 2) {
        // Pre-shift: Point=(1,0) → Index=1
        // Post-shift: Point=(0,0) → Index=0
        densityField.translate({1, 2});

        REQUIRE_THAT(densityField.getDensity(0), Catch::Matchers::WithinAbs(0.4, NUMERICAL_ZERO));
    }

    else {
        // Pre-shift: Point=(1,0,0) → Index=1
        // Post-shift: Point=(0,0,1) → Index=4
        densityField.translate({1, 2, 3});

        REQUIRE_THAT(densityField.getDensity(4), Catch::Matchers::WithinAbs(0.4, NUMERICAL_ZERO));
    }
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test operator==", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField1(resolution);
    densityField1.setConstant(0.4);

    TestType densityField2(resolution);
    densityField2.setConstant(0.4);

    REQUIRE(densityField1 == densityField2);
}

TEMPLATE_LIST_TEST_CASE("monad::DensityField: Test operator!=", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(2);

    TestType densityField1(resolution);
    densityField1.setConstant(0.4);

    TestType densityField2(resolution);
    densityField2.setConstant(0.5);

    REQUIRE(densityField1 != densityField2);
}
