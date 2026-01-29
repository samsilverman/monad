#include <tuple>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/detail/constants.hpp"
#include "test_common.hpp"

using namespace monad;
using namespace monad::testing;

using Types = std::tuple<Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test initalization", "[monad]", Types) {
    SECTION("Invalid resolution") {
        REQUIRE_THROWS_AS(TestType({0, 3, 4}, {0.5, 1.5, 2.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 0, 4}, {0.5, 1.5, 2.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3, 0}, {0.5, 1.5, 2.0}), std::invalid_argument);
    }

    SECTION("Invalid size") {
        REQUIRE_THROWS_AS(TestType({2, 3, 4}, {0.0, 1.5, 2.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3, 4}, {-0.5, 1.5, 2.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3, 4}, {0.5, 0.0, 2.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3, 4}, {0.5, -1.5, 2.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3, 4}, {0.5, 1.5, 0.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3, 4}, {0.5, 1.5, -2.0}), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test resolution", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    const TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    const Resolution expected{2, 3, 4};

    REQUIRE(grid.resolution() == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test size", "[monad]", Types) {
    using Size = typename TestType::Size;

    const TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    const Size expected{0.5, 1.5, 2.0};

    REQUIRE(grid.size() == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test densities", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    const TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    const DensityList expected(24, 0.0);

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test numElements", "[monad]", Types) {
    const TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    REQUIRE(grid.numElements() == 24);
}

TEST_CASE("monad::Grid3d: Test numNodes", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        REQUIRE(grid.numNodes() == 60);
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        REQUIRE(grid.numNodes() == 193);
    }
}

TEST_CASE("monad::Grid3d: Test numPeriodicNodes", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        REQUIRE(grid.numPeriodicNodes() == 24);
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        REQUIRE(grid.numPeriodicNodes() == 96);
    }
}

TEST_CASE("monad::Grid3d: Test node", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex8Grid::Point expected(0.25, 0.0, 0.0);

            REQUIRE(grid.node(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.node(grid.numNodes()), std::out_of_range);
        }
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex20Grid::Point expected(0.25, 0.0, 0.0);

            REQUIRE(grid.node(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.node(grid.numNodes()), std::out_of_range);
        }
    }
}

TEST_CASE("monad::Grid3d: Test nodes", "[monad]") {
    SECTION("Hex8Grid") {
        using NodesList = typename Hex8Grid::NodesList;

        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        const auto actual = grid.nodes();

        const NodesList expected {
            {0.0, 0.0, 0.0},
            {0.25, 0.0, 0.0},
            {0.5, 0.0, 0.0},
            {0.0, 0.5, 0.0},
            {0.25, 0.5, 0.0},
            {0.5, 0.5, 0.0},
            {0.0, 1.0, 0.0},
            {0.25, 1.0, 0.0},
            {0.5, 1.0, 0.0},
            {0.0, 1.5, 0.0},
            {0.25, 1.5, 0.0},
            {0.5, 1.5, 0.0},
            {0.0, 0.0, 0.5},
            {0.25, 0.0, 0.5},
            {0.5, 0.0, 0.5},
            {0.0, 0.5, 0.5},
            {0.25, 0.5, 0.5},
            {0.5, 0.5, 0.5},
            {0.0, 1.0, 0.5},
            {0.25, 1.0, 0.5},
            {0.5, 1.0, 0.5},
            {0.0, 1.5, 0.5},
            {0.25, 1.5, 0.5},
            {0.5, 1.5, 0.5},
            {0.0, 0.0, 1.0},
            {0.25, 0.0, 1.0},
            {0.5, 0.0, 1.0},
            {0.0, 0.5, 1.0},
            {0.25, 0.5, 1.0},
            {0.5, 0.5, 1.0},
            {0.0, 1.0, 1.0},
            {0.25, 1.0, 1.0},
            {0.5, 1.0, 1.0},
            {0.0, 1.5, 1.0},
            {0.25, 1.5, 1.0},
            {0.5, 1.5, 1.0},
            {0.0, 0.0, 1.5},
            {0.25, 0.0, 1.5},
            {0.5, 0.0, 1.5},
            {0.0, 0.5, 1.5},
            {0.25, 0.5, 1.5},
            {0.5, 0.5, 1.5},
            {0.0, 1.0, 1.5},
            {0.25, 1.0, 1.5},
            {0.5, 1.0, 1.5},
            {0.0, 1.5, 1.5},
            {0.25, 1.5, 1.5},
            {0.5, 1.5, 1.5},
            {0.0, 0.0, 2.0},
            {0.25, 0.0, 2.0},
            {0.5, 0.0, 2.0},
            {0.0, 0.5, 2.0},
            {0.25, 0.5, 2.0},
            {0.5, 0.5, 2.0},
            {0.0, 1.0, 2.0},
            {0.25, 1.0, 2.0},
            {0.5, 1.0, 2.0},
            {0.0, 1.5, 2.0},
            {0.25, 1.5, 2.0},
            {0.5, 1.5, 2.0}
        };

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }

    SECTION("Hex20Grid") {
        using NodesList = typename Hex20Grid::NodesList;

        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        const auto actual = grid.nodes();

        const NodesList expected {
            {0.0, 0.0, 0.0},
            {0.25, 0.0, 0.0},
            {0.5, 0.0, 0.0},
            {0.0, 0.5, 0.0},
            {0.25, 0.5, 0.0},
            {0.5, 0.5, 0.0},
            {0.0, 1.0, 0.0},
            {0.25, 1.0, 0.0},
            {0.5, 1.0, 0.0},
            {0.0, 1.5, 0.0},
            {0.25, 1.5, 0.0},
            {0.5, 1.5, 0.0},
            {0.0, 0.0, 0.5},
            {0.25, 0.0, 0.5},
            {0.5, 0.0, 0.5},
            {0.0, 0.5, 0.5},
            {0.25, 0.5, 0.5},
            {0.5, 0.5, 0.5},
            {0.0, 1.0, 0.5},
            {0.25, 1.0, 0.5},
            {0.5, 1.0, 0.5},
            {0.0, 1.5, 0.5},
            {0.25, 1.5, 0.5},
            {0.5, 1.5, 0.5},
            {0.0, 0.0, 1.0},
            {0.25, 0.0, 1.0},
            {0.5, 0.0, 1.0},
            {0.0, 0.5, 1.0},
            {0.25, 0.5, 1.0},
            {0.5, 0.5, 1.0},
            {0.0, 1.0, 1.0},
            {0.25, 1.0, 1.0},
            {0.5, 1.0, 1.0},
            {0.0, 1.5, 1.0},
            {0.25, 1.5, 1.0},
            {0.5, 1.5, 1.0},
            {0.0, 0.0, 1.5},
            {0.25, 0.0, 1.5},
            {0.5, 0.0, 1.5},
            {0.0, 0.5, 1.5},
            {0.25, 0.5, 1.5},
            {0.5, 0.5, 1.5},
            {0.0, 1.0, 1.5},
            {0.25, 1.0, 1.5},
            {0.5, 1.0, 1.5},
            {0.0, 1.5, 1.5},
            {0.25, 1.5, 1.5},
            {0.5, 1.5, 1.5},
            {0.0, 0.0, 2.0},
            {0.25, 0.0, 2.0},
            {0.5, 0.0, 2.0},
            {0.0, 0.5, 2.0},
            {0.25, 0.5, 2.0},
            {0.5, 0.5, 2.0},
            {0.0, 1.0, 2.0},
            {0.25, 1.0, 2.0},
            {0.5, 1.0, 2.0},
            {0.0, 1.5, 2.0},
            {0.25, 1.5, 2.0},
            {0.5, 1.5, 2.0},
            {0.125, 0.0, 0.0},
            {0.375, 0.0, 0.0},
            {0.125, 0.5, 0.0},
            {0.375, 0.5, 0.0},
            {0.125, 1.0, 0.0},
            {0.375, 1.0, 0.0},
            {0.125, 1.5, 0.0},
            {0.375, 1.5, 0.0},
            {0.125, 0.0, 0.5},
            {0.375, 0.0, 0.5},
            {0.125, 0.5, 0.5},
            {0.375, 0.5, 0.5},
            {0.125, 1.0, 0.5},
            {0.375, 1.0, 0.5},
            {0.125, 1.5, 0.5},
            {0.375, 1.5, 0.5},
            {0.125, 0.0, 1.0},
            {0.375, 0.0, 1.0},
            {0.125, 0.5, 1.0},
            {0.375, 0.5, 1.0},
            {0.125, 1.0, 1.0},
            {0.375, 1.0, 1.0},
            {0.125, 1.5, 1.0},
            {0.375, 1.5, 1.0},
            {0.125, 0.0, 1.5},
            {0.375, 0.0, 1.5},
            {0.125, 0.5, 1.5},
            {0.375, 0.5, 1.5},
            {0.125, 1.0, 1.5},
            {0.375, 1.0, 1.5},
            {0.125, 1.5, 1.5},
            {0.375, 1.5, 1.5},
            {0.125, 0.0, 2.0},
            {0.375, 0.0, 2.0},
            {0.125, 0.5, 2.0},
            {0.375, 0.5, 2.0},
            {0.125, 1.0, 2.0},
            {0.375, 1.0, 2.0},
            {0.125, 1.5, 2.0},
            {0.375, 1.5, 2.0},
            {0.0, 0.25, 0.0},
            {0.25, 0.25, 0.0},
            {0.5, 0.25, 0.0},
            {0.0, 0.75, 0.0},
            {0.25, 0.75, 0.0},
            {0.5, 0.75, 0.0},
            {0.0, 1.25, 0.0},
            {0.25, 1.25, 0.0},
            {0.5, 1.25, 0.0},
            {0.0, 0.25, 0.5},
            {0.25, 0.25, 0.5},
            {0.5, 0.25, 0.5},
            {0.0, 0.75, 0.5},
            {0.25, 0.75, 0.5},
            {0.5, 0.75, 0.5},
            {0.0, 1.25, 0.5},
            {0.25, 1.25, 0.5},
            {0.5, 1.25, 0.5},
            {0.0, 0.25, 1.0},
            {0.25, 0.25, 1.0},
            {0.5, 0.25, 1.0},
            {0.0, 0.75, 1.0},
            {0.25, 0.75, 1.0},
            {0.5, 0.75, 1.0},
            {0.0, 1.25, 1.0},
            {0.25, 1.25, 1.0},
            {0.5, 1.25, 1.0},
            {0.0, 0.25, 1.5},
            {0.25, 0.25, 1.5},
            {0.5, 0.25, 1.5},
            {0.0, 0.75, 1.5},
            {0.25, 0.75, 1.5},
            {0.5, 0.75, 1.5},
            {0.0, 1.25, 1.5},
            {0.25, 1.25, 1.5},
            {0.5, 1.25, 1.5},
            {0.0, 0.25, 2.0},
            {0.25, 0.25, 2.0},
            {0.5, 0.25, 2.0},
            {0.0, 0.75, 2.0},
            {0.25, 0.75, 2.0},
            {0.5, 0.75, 2.0},
            {0.0, 1.25, 2.0},
            {0.25, 1.25, 2.0},
            {0.5, 1.25, 2.0},
            {0.0, 0.0, 0.25},
            {0.25, 0.0, 0.25},
            {0.5, 0.0, 0.25},
            {0.0, 0.5, 0.25},
            {0.25, 0.5, 0.25},
            {0.5, 0.5, 0.25},
            {0.0, 1.0, 0.25},
            {0.25, 1.0, 0.25},
            {0.5, 1.0, 0.25},
            {0.0, 1.5, 0.25},
            {0.25, 1.5, 0.25},
            {0.5, 1.5, 0.25},
            {0.0, 0.0, 0.75},
            {0.25, 0.0, 0.75},
            {0.5, 0.0, 0.75},
            {0.0, 0.5, 0.75},
            {0.25, 0.5, 0.75},
            {0.5, 0.5, 0.75},
            {0.0, 1.0, 0.75},
            {0.25, 1.0, 0.75},
            {0.5, 1.0, 0.75},
            {0.0, 1.5, 0.75},
            {0.25, 1.5, 0.75},
            {0.5, 1.5, 0.75},
            {0.0, 0.0, 1.25},
            {0.25, 0.0, 1.25},
            {0.5, 0.0, 1.25},
            {0.0, 0.5, 1.25},
            {0.25, 0.5, 1.25},
            {0.5, 0.5, 1.25},
            {0.0, 1.0, 1.25},
            {0.25, 1.0, 1.25},
            {0.5, 1.0, 1.25},
            {0.0, 1.5, 1.25},
            {0.25, 1.5, 1.25},
            {0.5, 1.5, 1.25},
            {0.0, 0.0, 1.75},
            {0.25, 0.0, 1.75},
            {0.5, 0.0, 1.75},
            {0.0, 0.5, 1.75},
            {0.25, 0.5, 1.75},
            {0.5, 0.5, 1.75},
            {0.0, 1.0, 1.75},
            {0.25, 1.0, 1.75},
            {0.5, 1.0, 1.75},
            {0.0, 1.5, 1.75},
            {0.25, 1.5, 1.75},
            {0.5, 1.5, 1.75}
        };

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }
}

TEST_CASE("monad::Grid3d: Test element", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex8Grid::ElementList expected{1, 2, 5, 4, 13, 14, 17, 16};

            REQUIRE(grid.element(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.element(grid.numElements()), std::out_of_range);
        }
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex20Grid::ElementList expected{1, 2, 5, 4, 13, 14, 17, 16, 61, 102, 63, 101, 69, 111, 71, 110, 146, 147, 150, 149};

            REQUIRE(grid.element(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.element(grid.numElements()), std::out_of_range);
        }
    }
}

TEST_CASE("monad::Grid3d: Test elements", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        const Hex8Grid::ElementsList expected {
            {0, 1, 4, 3, 12, 13, 16, 15},
            {1, 2, 5, 4, 13, 14, 17, 16},
            {3, 4, 7, 6, 15, 16, 19, 18},
            {4, 5, 8, 7, 16, 17, 20, 19},
            {6, 7, 10, 9, 18, 19, 22, 21},
            {7, 8, 11, 10, 19, 20, 23, 22},
            {12, 13, 16, 15, 24, 25, 28, 27},
            {13, 14, 17, 16, 25, 26, 29, 28},
            {15, 16, 19, 18, 27, 28, 31, 30},
            {16, 17, 20, 19, 28, 29, 32, 31},
            {18, 19, 22, 21, 30, 31, 34, 33},
            {19, 20, 23, 22, 31, 32, 35, 34},
            {24, 25, 28, 27, 36, 37, 40, 39},
            {25, 26, 29, 28, 37, 38, 41, 40},
            {27, 28, 31, 30, 39, 40, 43, 42},
            {28, 29, 32, 31, 40, 41, 44, 43},
            {30, 31, 34, 33, 42, 43, 46, 45},
            {31, 32, 35, 34, 43, 44, 47, 46},
            {36, 37, 40, 39, 48, 49, 52, 51},
            {37, 38, 41, 40, 49, 50, 53, 52},
            {39, 40, 43, 42, 51, 52, 55, 54},
            {40, 41, 44, 43, 52, 53, 56, 55},
            {42, 43, 46, 45, 54, 55, 58, 57},
            {43, 44, 47, 46, 55, 56, 59, 58}
        };

        REQUIRE(grid.elements() == expected);
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        const Hex20Grid::ElementsList expected {
            {0, 1, 4, 3, 12, 13, 16, 15, 60, 101, 62, 100, 68, 110, 70, 109, 145, 146, 149, 148},
            {1, 2, 5, 4, 13, 14, 17, 16, 61, 102, 63, 101, 69, 111, 71, 110, 146, 147, 150, 149},
            {3, 4, 7, 6, 15, 16, 19, 18, 62, 104, 64, 103, 70, 113, 72, 112, 148, 149, 152, 151},
            {4, 5, 8, 7, 16, 17, 20, 19, 63, 105, 65, 104, 71, 114, 73, 113, 149, 150, 153, 152},
            {6, 7, 10, 9, 18, 19, 22, 21, 64, 107, 66, 106, 72, 116, 74, 115, 151, 152, 155, 154},
            {7, 8, 11, 10, 19, 20, 23, 22, 65, 108, 67, 107, 73, 117, 75, 116, 152, 153, 156, 155},
            {12, 13, 16, 15, 24, 25, 28, 27, 68, 110, 70, 109, 76, 119, 78, 118, 157, 158, 161, 160},
            {13, 14, 17, 16, 25, 26, 29, 28, 69, 111, 71, 110, 77, 120, 79, 119, 158, 159, 162, 161},
            {15, 16, 19, 18, 27, 28, 31, 30, 70, 113, 72, 112, 78, 122, 80, 121, 160, 161, 164, 163},
            {16, 17, 20, 19, 28, 29, 32, 31, 71, 114, 73, 113, 79, 123, 81, 122, 161, 162, 165, 164},
            {18, 19, 22, 21, 30, 31, 34, 33, 72, 116, 74, 115, 80, 125, 82, 124, 163, 164, 167, 166},
            {19, 20, 23, 22, 31, 32, 35, 34, 73, 117, 75, 116, 81, 126, 83, 125, 164, 165, 168, 167},
            {24, 25, 28, 27, 36, 37, 40, 39, 76, 119, 78, 118, 84, 128, 86, 127, 169, 170, 173, 172},
            {25, 26, 29, 28, 37, 38, 41, 40, 77, 120, 79, 119, 85, 129, 87, 128, 170, 171, 174, 173},
            {27, 28, 31, 30, 39, 40, 43, 42, 78, 122, 80, 121, 86, 131, 88, 130, 172, 173, 176, 175},
            {28, 29, 32, 31, 40, 41, 44, 43, 79, 123, 81, 122, 87, 132, 89, 131, 173, 174, 177, 176},
            {30, 31, 34, 33, 42, 43, 46, 45, 80, 125, 82, 124, 88, 134, 90, 133, 175, 176, 179, 178},
            {31, 32, 35, 34, 43, 44, 47, 46, 81, 126, 83, 125, 89, 135, 91, 134, 176, 177, 180, 179},
            {36, 37, 40, 39, 48, 49, 52, 51, 84, 128, 86, 127, 92, 137, 94, 136, 181, 182, 185, 184},
            {37, 38, 41, 40, 49, 50, 53, 52, 85, 129, 87, 128, 93, 138, 95, 137, 182, 183, 186, 185},
            {39, 40, 43, 42, 51, 52, 55, 54, 86, 131, 88, 130, 94, 140, 96, 139, 184, 185, 188, 187},
            {40, 41, 44, 43, 52, 53, 56, 55, 87, 132, 89, 131, 95, 141, 97, 140, 185, 186, 189, 188},
            {42, 43, 46, 45, 54, 55, 58, 57, 88, 134, 90, 133, 96, 143, 98, 142, 187, 188, 191, 190},
            {43, 44, 47, 46, 55, 56, 59, 58, 89, 135, 91, 134, 97, 144, 99, 143, 188, 189, 192, 191}
        };

        REQUIRE(grid.elements() == expected);
    }
}

TEST_CASE("monad::Grid3d: Test periodicElement", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex8Grid::ElementList expected{1, 0, 2, 3, 7, 6, 8, 9};

            REQUIRE(grid.periodicElement(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.periodicElement(grid.numElements()), std::out_of_range);
        }
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex20Grid::ElementList expected{1, 0, 2, 3, 7, 6, 8, 9, 25, 48, 27, 49, 31, 54, 33, 55, 73, 72, 74, 75};

            REQUIRE(grid.periodicElement(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.periodicElement(grid.numElements()), std::out_of_range);
        }
    }
}

TEST_CASE("monad::Grid3d: Test periodicElements", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        const Hex8Grid::ElementsList expected {
            {0, 1, 3, 2, 6, 7, 9, 8},
            {1, 0, 2, 3, 7, 6, 8, 9},
            {2, 3, 5, 4, 8, 9, 11, 10},
            {3, 2, 4, 5, 9, 8, 10, 11},
            {4, 5, 1, 0, 10, 11, 7, 6},
            {5, 4, 0, 1, 11, 10, 6, 7},
            {6, 7, 9, 8, 12, 13, 15, 14},
            {7, 6, 8, 9, 13, 12, 14, 15},
            {8, 9, 11, 10, 14, 15, 17, 16},
            {9, 8, 10, 11, 15, 14, 16, 17},
            {10, 11, 7, 6, 16, 17, 13, 12},
            {11, 10, 6, 7, 17, 16, 12, 13},
            {12, 13, 15, 14, 18, 19, 21, 20},
            {13, 12, 14, 15, 19, 18, 20, 21},
            {14, 15, 17, 16, 20, 21, 23, 22},
            {15, 14, 16, 17, 21, 20, 22, 23},
            {16, 17, 13, 12, 22, 23, 19, 18},
            {17, 16, 12, 13, 23, 22, 18, 19},
            {18, 19, 21, 20, 0, 1, 3, 2},
            {19, 18, 20, 21, 1, 0, 2, 3},
            {20, 21, 23, 22, 2, 3, 5, 4},
            {21, 20, 22, 23, 3, 2, 4, 5},
            {22, 23, 19, 18, 4, 5, 1, 0},
            {23, 22, 18, 19, 5, 4, 0, 1}
        };

        REQUIRE(grid.periodicElements() == expected);
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        const Hex20Grid::ElementsList expected {
            {0, 1, 3, 2, 6, 7, 9, 8, 24, 49, 26, 48, 30, 55, 32, 54, 72, 73, 75, 74},
            {1, 0, 2, 3, 7, 6, 8, 9, 25, 48, 27, 49, 31, 54, 33, 55, 73, 72, 74, 75},
            {2, 3, 5, 4, 8, 9, 11, 10, 26, 51, 28, 50, 32, 57, 34, 56, 74, 75, 77, 76},
            {3, 2, 4, 5, 9, 8, 10, 11, 27, 50, 29, 51, 33, 56, 35, 57, 75, 74, 76, 77},
            {4, 5, 1, 0, 10, 11, 7, 6, 28, 53, 24, 52, 34, 59, 30, 58, 76, 77, 73, 72},
            {5, 4, 0, 1, 11, 10, 6, 7, 29, 52, 25, 53, 35, 58, 31, 59, 77, 76, 72, 73},
            {6, 7, 9, 8, 12, 13, 15, 14, 30, 55, 32, 54, 36, 61, 38, 60, 78, 79, 81, 80},
            {7, 6, 8, 9, 13, 12, 14, 15, 31, 54, 33, 55, 37, 60, 39, 61, 79, 78, 80, 81},
            {8, 9, 11, 10, 14, 15, 17, 16, 32, 57, 34, 56, 38, 63, 40, 62, 80, 81, 83, 82},
            {9, 8, 10, 11, 15, 14, 16, 17, 33, 56, 35, 57, 39, 62, 41, 63, 81, 80, 82, 83},
            {10, 11, 7, 6, 16, 17, 13, 12, 34, 59, 30, 58, 40, 65, 36, 64, 82, 83, 79, 78},
            {11, 10, 6, 7, 17, 16, 12, 13, 35, 58, 31, 59, 41, 64, 37, 65, 83, 82, 78, 79},
            {12, 13, 15, 14, 18, 19, 21, 20, 36, 61, 38, 60, 42, 67, 44, 66, 84, 85, 87, 86},
            {13, 12, 14, 15, 19, 18, 20, 21, 37, 60, 39, 61, 43, 66, 45, 67, 85, 84, 86, 87},
            {14, 15, 17, 16, 20, 21, 23, 22, 38, 63, 40, 62, 44, 69, 46, 68, 86, 87, 89, 88},
            {15, 14, 16, 17, 21, 20, 22, 23, 39, 62, 41, 63, 45, 68, 47, 69, 87, 86, 88, 89},
            {16, 17, 13, 12, 22, 23, 19, 18, 40, 65, 36, 64, 46, 71, 42, 70, 88, 89, 85, 84},
            {17, 16, 12, 13, 23, 22, 18, 19, 41, 64, 37, 65, 47, 70, 43, 71, 89, 88, 84, 85},
            {18, 19, 21, 20, 0, 1, 3, 2, 42, 67, 44, 66, 24, 49, 26, 48, 90, 91, 93, 92},
            {19, 18, 20, 21, 1, 0, 2, 3, 43, 66, 45, 67, 25, 48, 27, 49, 91, 90, 92, 93},
            {20, 21, 23, 22, 2, 3, 5, 4, 44, 69, 46, 68, 26, 51, 28, 50, 92, 93, 95, 94},
            {21, 20, 22, 23, 3, 2, 4, 5, 45, 68, 47, 69, 27, 50, 29, 51, 93, 92, 94, 95},
            {22, 23, 19, 18, 4, 5, 1, 0, 46, 71, 42, 70, 28, 53, 24, 52, 94, 95, 91, 90},
            {23, 22, 18, 19, 5, 4, 0, 1, 47, 70, 43, 71, 29, 52, 25, 53, 95, 94, 90, 91}
        };

        REQUIRE(grid.periodicElements() == expected);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test getDensity/setDensity", "[monad]", Types) {
    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    SECTION("No errors") {
        REQUIRE_THAT(grid.getDensity(1), Catch::Matchers::WithinAbs(0.0, NUMERICAL_ZERO));

        grid.setDensity(1, 0.1);

        REQUIRE_THAT(grid.getDensity(1), Catch::Matchers::WithinAbs(0.1, NUMERICAL_ZERO));
    }

    SECTION("Invalid index") {
        REQUIRE_THROWS_AS(grid.getDensity(grid.numElements()), std::out_of_range);

        REQUIRE_THROWS_AS(grid.setDensity(grid.numElements(), 0.1), std::out_of_range);
    }

    SECTION("Invalid density") {
        REQUIRE_THROWS_AS(grid.setDensity(1, -0.1), std::invalid_argument);
        REQUIRE_THROWS_AS(grid.setDensity(1, 1.1), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test setDensities", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    SECTION("No errors") {
        const DensityList densities(grid.numElements(), 0.5);

        grid.setDensities(densities);

        REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(densities).margin(NUMERICAL_ZERO));
    }

    SECTION("Invalid size") {
        const DensityList densities(grid.numElements() + 1, 0.5);

        REQUIRE_THROWS_AS(grid.setDensities(densities), std::invalid_argument);
    }

    SECTION("Invalid density") {
        DensityList densities(grid.numElements(), 0.5);
        densities[1] *= -1;

        REQUIRE_THROWS_AS(grid.setDensities(densities), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test setDensitiesConstant", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    SECTION("No errors") {
        grid.setDensitiesConstant(0.5);

        const DensityList expected(grid.numElements(), 0.5);

        REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
    }

    SECTION("Invalid density") {
        REQUIRE_THROWS_AS(grid.setDensitiesConstant(-0.1), std::invalid_argument);
        REQUIRE_THROWS_AS(grid.setDensitiesConstant(1.1), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test setDensitiesZeros", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    grid.setDensitiesZeros();

    const DensityList expected(grid.numElements(), 0.0);

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test setDensitiesOnes", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    grid.setDensitiesOnes();

    const DensityList expected(grid.numElements(), 1.0);

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test setDensitiesRandom", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    grid.setDensitiesRandom(1234);

    const DensityList expected {
        0.49766366630595177,
        0.81783844288953345,
        0.61211189358442286,
        0.77135991905475021,
        0.86066977348643747,
        0.15063696570521984,
        0.19851876010104652,
        0.8151629340839277,
        0.15881535340057207,
        0.11613783027882053,
        0.01290753305838846,
        0.48683344491314973,
        0.33101542682467566,
        0.80263957374559314,
        0.09825193844042972,
        0.05599344953922467,
        0.44266265526020054,
        0.02214390901237923,
        0.29072854825594158,
        0.24639444143572856,
        0.73828676607573551,
        0.88922613337994483,
        0.98713930325475108,
        0.1174433832915556
    };

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test setDensitiesFunction", "[monad]", Types) {
    using Point = typename TestType::Point;

    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    SECTION("No errors") {
        // sin²(x+y+z)
        auto f = [](const Point &point) -> double {
            const double x = point(0);
            const double y = point(1);
            const double z = point(2);

            return std::pow(std::sin(x + y + z), 2);
        };

        grid.setDensitiesFunction(f);

        double expected = 0.0;
        const auto nodes = grid.elementNodes(1);

        for (Point node : nodes.rowwise()) {
            expected += f(node);
        }

        expected /= static_cast<double>(nodes.rows());

        REQUIRE_THAT(grid.getDensity(1), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
    }

    SECTION("Bad function") {
        // exp(x+y+z)
        auto f = [](const Point &point) -> double {
            const double x = point(0);
            const double y = point(1);
            const double z = point(2);

            return std::exp(x + y + z);
        };

        REQUIRE_THROWS_AS(grid.setDensitiesFunction(f), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test translate", "[monad]", Types) {
    TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    // Index=1 → Point=(1,0,0)
    grid.setDensity(1, 0.5);

    // Point=(0,2,3)
    grid.translate({1, 2, 3});

    REQUIRE_THAT(grid.getDensity(22), Catch::Matchers::WithinAbs(0.5, NUMERICAL_ZERO));
}

TEST_CASE("monad::Grid3d: Test elementNodes", "[monad]") {
    SECTION("Hex8Grid") {
        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex8Grid::NodesMatrix expected {
                {0.25, 0.0, 0.0},
                {0.5, 0.0, 0.0},
                {0.5, 0.5, 0.0},
                {0.25, 0.5, 0.0},
                {0.25, 0.0, 0.5},
                {0.5, 0.0, 0.5},
                {0.5, 0.5, 0.5},
                {0.25, 0.5, 0.5}
            };

            REQUIRE(grid.elementNodes(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.elementNodes(grid.numElements()), std::out_of_range);
        }
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

        SECTION("No errors") {
            const Hex20Grid::NodesMatrix expected {
                {0.25, 0.0, 0.0},
                {0.5, 0.0, 0.0},
                {0.5, 0.5, 0.0},
                {0.25, 0.5, 0.0},
                {0.25, 0.0, 0.5},
                {0.5, 0.0, 0.5},
                {0.5, 0.5, 0.5},
                {0.25, 0.5, 0.5},
                {0.375, 0.0, 0.0},
                {0.5, 0.25, 0.0},
                {0.375, 0.5, 0.0},
                {0.25, 0.25, 0.0},
                {0.375, 0.0, 0.5},
                {0.5, 0.25, 0.5},
                {0.375, 0.5, 0.5},
                {0.25, 0.25, 0.5},
                {0.25, 0.0, 0.25},
                {0.5, 0.0, 0.25},
                {0.5, 0.5, 0.25},
                {0.25, 0.5, 0.25}
            };

            REQUIRE(grid.elementNodes(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.elementNodes(grid.numElements()), std::out_of_range);
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test measure/volume", "[monad]", Types) {
    const TestType grid({2, 3, 4}, {0.5, 1.5, 2.0});

    REQUIRE_THAT(grid.measure(), Catch::Matchers::WithinAbs(1.5, NUMERICAL_ZERO));
    REQUIRE_THAT(grid.volume(), Catch::Matchers::WithinAbs(1.5, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test operator==", "[monad]", Types) {
    const TestType grid1({2, 3, 4}, {0.5, 1.5, 2.0});
    const TestType grid2({2, 3, 4}, {0.5, 1.5, 2.0});

    REQUIRE(grid1 == grid2);
}

TEMPLATE_LIST_TEST_CASE("monad::Grid3d: Test operator!=", "[monad]", Types) {
    const TestType grid1({2, 3, 4}, {0.5, 1.5, 2.0});

    SECTION("Different resolution") {
        const TestType grid2({3, 3, 4}, {0.5, 1.5, 2.0});
        const TestType grid3({2, 4, 4}, {0.5, 1.5, 2.0});
        const TestType grid4({2, 3, 5}, {0.5, 1.5, 2.0});

        REQUIRE(grid1 != grid2);
        REQUIRE(grid1 != grid3);
        REQUIRE(grid1 != grid4);
    }

    SECTION("Different size") {
        const TestType grid2({2, 3, 4}, {0.6, 1.5, 2.0});
        const TestType grid3({2, 3, 4}, {0.5, 1.6, 2.0});
        const TestType grid4({2, 3, 4}, {0.5, 1.5, 2.1});

        REQUIRE(grid1 != grid2);
        REQUIRE(grid1 != grid3);
        REQUIRE(grid1 != grid4);
    }
}
