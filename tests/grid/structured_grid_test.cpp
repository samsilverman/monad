#include <tuple>
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test initalization", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    SECTION("Invalid resolution") {
        for (std::size_t i = 0; i < TestType::Dim; ++i) {
            Resolution resolutionInvalid = resolution;
            resolutionInvalid[i] = 0;

            REQUIRE_THROWS_AS(TestType(resolutionInvalid, size), std::invalid_argument);
        }
    }

    SECTION("Invalid size") {
        for (std::size_t i = 0; i < TestType::Dim; ++i) {
            Size sizeInvalid = size;
            sizeInvalid[i] = 0.0;

            REQUIRE_THROWS_AS(TestType(resolution, sizeInvalid), std::invalid_argument);

            sizeInvalid[i] = -1.0;

            REQUIRE_THROWS_AS(TestType(resolution, sizeInvalid), std::invalid_argument);
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test resolution", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    REQUIRE(grid.resolution() == resolution);
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test size", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    REQUIRE(grid.size() == size);
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test numElements", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    const std::size_t expected = static_cast<std::size_t>(std::pow(2.0, TestType::Dim));
    const std::size_t actual = grid.numElements();

    REQUIRE(actual == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test numNodes", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Topology = typename TestType::Topology;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    const std::size_t expected = Topology::numNodes(resolution);
    const std::size_t actual = grid.numNodes();

    REQUIRE(actual == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test numPeriodicNodes", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Topology = typename TestType::Topology;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    const std::size_t expected = Topology::numPeriodicNodes(resolution);
    const std::size_t actual = grid.numPeriodicNodes();

    REQUIRE(actual == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test node", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Topology = typename TestType::Topology;
    using Point = typename TestType::Point;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    SECTION("No errors") {
        const Point expected = Topology::node(1, resolution, size);
        const Point actual = grid.node(1);

        REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
    }

    SECTION("Invalid index") {
        const std::size_t badIndex = grid.numNodes();
    
        REQUIRE_THROWS_AS(grid.node(badIndex), std::out_of_range);
    }
}

TEST_CASE("monad::grid::StructuredGrid: Test nodes", "[monad]") {
    SECTION("Quad4Grid") {
        using NodesList = typename Quad4Grid::NodesList;

        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        const NodesList expected {
            {0.0, 0.0},
            {0.25, 0.0},
            {0.5, 0.0},
            {0.0, 0.5},
            {0.25, 0.5},
            {0.5, 0.5},
            {0.0, 1.0},
            {0.25, 1.0},
            {0.5, 1.0},
            {0.0, 1.5},
            {0.25, 1.5},
            {0.5, 1.5}
        };

        const auto actual = grid.nodes();

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }

    SECTION("Quad8Grid") {
        using NodesList = typename Quad8Grid::NodesList;

        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        const NodesList expected {
            {0.0, 0.0},
            {0.25, 0.0},
            {0.5, 0.0},
            {0.0, 0.5},
            {0.25, 0.5},
            {0.5, 0.5},
            {0.0, 1.0},
            {0.25, 1.0},
            {0.5, 1.0},
            {0.0, 1.5},
            {0.25, 1.5},
            {0.5, 1.5},
            {0.125, 0.0},
            {0.375, 0.0},
            {0.125, 0.5},
            {0.375, 0.5},
            {0.125, 1.0},
            {0.375, 1.0},
            {0.125, 1.5},
            {0.375, 1.5},
            {0.0, 0.25},
            {0.25, 0.25},
            {0.5, 0.25},
            {0.0, 0.75},
            {0.25, 0.75},
            {0.5, 0.75},
            {0.0, 1.25},
            {0.25, 1.25},
            {0.5, 1.25}
        };

        const auto actual = grid.nodes();

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }

    SECTION("Hex8Grid") {
        using NodesList = typename Hex8Grid::NodesList;

        const Hex8Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

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

        const auto actual = grid.nodes();

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }

    SECTION("Hex20Grid") {
        using NodesList = typename Hex20Grid::NodesList;

        const Hex20Grid grid({2, 3, 4}, {0.5, 1.5, 2.0});

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

        const auto actual = grid.nodes();

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test element", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Topology = typename TestType::Topology;
    using ElementList = typename TestType::ElementList;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    SECTION("No errors") {
        const ElementList expected = Topology::element(1, resolution);
        const ElementList actual = grid.element(1);

        REQUIRE(actual == expected);
    }

    SECTION("Invalid index") {
        const std::size_t badIndex = grid.numElements();

        REQUIRE_THROWS_AS(grid.element(badIndex), std::out_of_range);
    }
}

TEST_CASE("monad::grid::StructuredGrid: Test elements", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        const Quad4Grid::ElementsList expected {
            {0, 1, 4, 3},
            {1, 2, 5, 4},
            {3, 4, 7, 6},
            {4, 5, 8, 7},
            {6, 7, 10, 9},
            {7, 8, 11, 10}
        };

        REQUIRE(grid.elements() == expected);
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        const Quad8Grid::ElementsList expected {
            {0, 1, 4, 3, 12, 21, 14, 20},
            {1, 2, 5, 4, 13, 22, 15, 21},
            {3, 4, 7, 6, 14, 24, 16, 23},
            {4, 5, 8, 7, 15, 25, 17, 24},
            {6, 7, 10, 9, 16, 27, 18, 26},
            {7, 8, 11, 10, 17, 28, 19, 27}
        };

        REQUIRE(grid.elements() == expected);
    }

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

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test periodicElement", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Topology = typename TestType::Topology;
    using ElementList = typename TestType::ElementList;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    SECTION("No errors") {
        const ElementList expected = Topology::periodicElement(1, resolution);
        const ElementList actual = grid.periodicElement(1);

        REQUIRE(actual == expected);
    }

    SECTION("Invalid index") {
        const std::size_t badIndex = grid.numElements();
        
        REQUIRE_THROWS_AS(grid.element(badIndex), std::out_of_range);
    }
}

TEST_CASE("monad::grid::StructuredGrid: Test periodicElements", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        const Quad4Grid::ElementsList expected {
            {0, 1, 3, 2},
            {1, 0, 2, 3},
            {2, 3, 5, 4},
            {3, 2, 4, 5},
            {4, 5, 1, 0},
            {5, 4, 0, 1}
        };

        REQUIRE(grid.periodicElements() == expected);
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        const Quad8Grid::ElementsList expected {
            {0, 1, 3, 2, 6, 13, 8, 12},
            {1, 0, 2, 3, 7, 12, 9, 13},
            {2, 3, 5, 4, 8, 15, 10, 14},
            {3, 2, 4, 5, 9, 14, 11, 15},
            {4, 5, 1, 0, 10, 17, 6, 16},
            {5, 4, 0, 1, 11, 16, 7, 17}
        };

        REQUIRE(grid.periodicElements() == expected);
    }

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

TEST_CASE("monad::grid::StructuredGrid: Test elementNodes", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad4Grid::NodesMatrix expected {
                {0.25, 0.0},
                {0.5, 0.0},
                {0.5, 0.5},
                {0.25, 0.5}
            };

            const auto actual = grid.elementNodes(1);

            REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            const std::size_t badIndex = grid.numElements();

            REQUIRE_THROWS_AS(grid.elementNodes(badIndex), std::out_of_range);
        }
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad8Grid::NodesMatrix expected {
                {0.25, 0.0},
                {0.5, 0.0},
                {0.5, 0.5},
                {0.25, 0.5},
                {0.375, 0.0},
                {0.5, 0.25},
                {0.375, 0.5},
                {0.25, 0.25}
            };

            const auto actual = grid.elementNodes(1);

            REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            const std::size_t badIndex = grid.numElements();

            REQUIRE_THROWS_AS(grid.elementNodes(badIndex), std::out_of_range);
        }
    }

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

            const auto actual = grid.elementNodes(1);

            REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            const std::size_t badIndex = grid.numElements();

            REQUIRE_THROWS_AS(grid.elementNodes(badIndex), std::out_of_range);
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

            const auto actual = grid.elementNodes(1);

            REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            const std::size_t badIndex = grid.numElements();

            REQUIRE_THROWS_AS(grid.elementNodes(badIndex), std::out_of_range);
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test measure/area/volume", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid(resolution, size);

    const double expected = std::pow(0.4, TestType::Dim);
    const double actual = grid.measure();

    REQUIRE_THAT(actual, Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));

    if constexpr (TestType::Dim == 2) {
        REQUIRE_THAT(grid.area(), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
    }
    else {
        REQUIRE_THAT(grid.volume(), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
    }
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test operator==", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid1(resolution, size);
    const TestType grid2(resolution, size);

    REQUIRE(grid1 == grid2);
}

TEMPLATE_LIST_TEST_CASE("monad::grid::StructuredGrid: Test operator!=", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.4);

    const TestType grid1(resolution, size);

    SECTION("Different resolution") {
        for (std::size_t i = 0; i < TestType::Dim; ++i) {
            Resolution resolutionDifferent = resolution;
            resolutionDifferent[i] += 1;

            const TestType grid2(resolutionDifferent, size);

            REQUIRE(grid1 != grid2);
        }
    }

    SECTION("Different size") {
        for (std::size_t i = 0; i < TestType::Dim; ++i) {
            Size sizeDifferent = size;
            sizeDifferent[i] += 0.1;

            const TestType grid2(resolution, sizeDifferent);

            REQUIRE(grid1 != grid2);
        }
    }
}
