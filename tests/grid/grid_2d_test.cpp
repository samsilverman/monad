#include <tuple>
#include <stdexcept>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/detail/constants.hpp"
#include "test_common.hpp"

using namespace monad;
using namespace monad::testing;

using Types = std::tuple<Quad4Grid, Quad8Grid>;

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test initalization", "[monad]", Types) {
    SECTION("Invalid resolution") {
        REQUIRE_THROWS_AS(TestType({0, 3}, {0.5, 1.5}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 0}, {0.5, 1.5}), std::invalid_argument);
    }

    SECTION("Invalid size") {
        REQUIRE_THROWS_AS(TestType({2, 3}, {0.0, 1.5}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3}, {-0.5, 1.5}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3}, {0.5, 0.0}), std::invalid_argument);
        REQUIRE_THROWS_AS(TestType({2, 3}, {0.5, -1.5}), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test resolution", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    const TestType grid({2, 3}, {0.5, 1.5});

    const Resolution expected{2, 3};

    REQUIRE(grid.resolution() == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test size", "[monad]", Types) {
    using Size = typename TestType::Size;

    const TestType grid({2, 3}, {0.5, 1.5});

    const Size expected{0.5, 1.5};

    REQUIRE(grid.size() == expected);
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test densities", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    const TestType grid({2, 3}, {0.5, 1.5});

    const DensityList expected(6, 0.0);

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test numElements", "[monad]", Types) {
    const TestType grid({2, 3}, {0.5, 1.5});

    REQUIRE(grid.numElements() == 6);
}

TEST_CASE("monad::Grid2d: Test numNodes", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        REQUIRE(grid.numNodes() == 12);
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        REQUIRE(grid.numNodes() == 29);
    }
}

TEST_CASE("monad::Grid2d: Test numPeriodicNodes", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        REQUIRE(grid.numPeriodicNodes() == 6);
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        REQUIRE(grid.numPeriodicNodes() == 18);
    }
}

TEST_CASE("monad::Grid2d: Test node", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad4Grid::Point expected(0.25, 0.0);

            REQUIRE(grid.node(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.node(grid.numNodes()), std::out_of_range);
        }
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad8Grid::Point expected(0.25, 0.0);

            REQUIRE(grid.node(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.node(grid.numNodes()), std::out_of_range);
        }
    }
}

TEST_CASE("monad::Grid2d: Test nodes", "[monad]") {
    SECTION("Quad4Grid") {
        using NodesList = typename Quad4Grid::NodesList;

        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        const auto actual = grid.nodes();

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

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }

    SECTION("Quad8Grid") {
        using NodesList = typename Quad8Grid::NodesList;

        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        const auto actual = grid.nodes();

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

        REQUIRE(actual.size() == expected.size());

        for (std::size_t i = 0; i < actual.size(); ++i) {
            REQUIRE(actual[i].isApprox(expected[i], NUMERICAL_ZERO));
        }
    }
}

TEST_CASE("monad::Grid2d: Test element", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad4Grid::ElementList expected{1, 2, 5, 4};

            REQUIRE(grid.element(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.element(grid.numElements()), std::out_of_range);
        }
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad8Grid::ElementList expected{1, 2, 5, 4, 13, 22, 15, 21};

            REQUIRE(grid.element(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.element(grid.numElements()), std::out_of_range);
        }
    }
}

TEST_CASE("monad::Grid2d: Test elements", "[monad]") {
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
}

TEST_CASE("monad::Grid2d: Test periodicElement", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad4Grid::ElementList expected{1, 0, 2, 3};

            REQUIRE(grid.periodicElement(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.periodicElement(grid.numElements()), std::out_of_range);
        }
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad8Grid::ElementList expected{1, 0, 2, 3, 7, 12, 9, 13};

            REQUIRE(grid.periodicElement(1) == expected);
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.periodicElement(grid.numElements()), std::out_of_range);
        }
    }
}

TEST_CASE("monad::Grid2d: Test periodicElements", "[monad]") {
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
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test getDensity/setDensity", "[monad]", Types) {
    TestType grid({2, 3}, {0.5, 1.5});

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

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensities", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3}, {0.5, 1.5});

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

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensitiesConstant", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3}, {0.5, 1.5});

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

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensitiesZeros", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3}, {0.5, 1.5});

    grid.setDensitiesZeros();

    const DensityList expected(grid.numElements(), 0.0);

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensitiesOnes", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3}, {0.5, 1.5});

    grid.setDensitiesOnes();

    const DensityList expected(grid.numElements(), 1.0);

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensitiesRandom", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3}, {0.5, 1.5});

    grid.setDensitiesRandom(1234);

    const DensityList expected {
        0.49766366630595177,
        0.81783844288953345,
        0.61211189358442286,
        0.77135991905475021,
        0.86066977348643747,
        0.15063696570521984
    };

    REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensitiesFunction", "[monad]", Types) {
    using Point = typename TestType::Point;

    TestType grid({2, 3}, {0.5, 1.5});

    SECTION("No errors") {
        // 0.1x+0.2y
        auto f = [](const Point &p) -> double {
            return 0.1 * p(0) + 0.2 * p(1);
        };

        grid.setDensitiesFunction(f);

        // The average of a linear function over a convex region equals its value at the centroid
        const auto nodes = grid.elementNodes(1);
        const Point centroid = nodes.colwise().mean();
        const double expected = f(centroid);

        REQUIRE_THAT(grid.getDensity(1), Catch::Matchers::WithinAbs(expected, NUMERICAL_ZERO));
    }

    SECTION("Bad function") {
        // exp(x+y)
        auto f = [](const Point &point) -> double {
            const double x = point(0);
            const double y = point(1);

            return std::exp(x + y);
        };

        REQUIRE_THROWS_AS(grid.setDensitiesFunction(f), std::invalid_argument);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test setDensitiesFile", "[monad]", Types) {
    using DensityList = typename TestType::DensityList;

    TestType grid({2, 3}, {0.5, 1.5});

    const auto folder = std::filesystem::path(__FILE__).parent_path() / "csv";

    SECTION("No errors") {
        const auto file = folder / "good.csv";

        grid.setDensitiesFile(file.string());

        const DensityList expected {
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6
        };

        REQUIRE_THAT(grid.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
    }

    SECTION("Invalid csv data") {
        // data < 0
        auto file = folder / "bad1.csv";

        REQUIRE_THROWS_AS(grid.setDensitiesFile(file.string()), std::runtime_error);

        // data > 1
        file = folder / "bad2.csv";

        REQUIRE_THROWS_AS(grid.setDensitiesFile(file.string()), std::runtime_error);

        // not enough data
        file = folder / "bad3.csv";

        REQUIRE_THROWS_AS(grid.setDensitiesFile(file.string()), std::runtime_error);

        // too much data
        file = folder / "bad4.csv";

        REQUIRE_THROWS_AS(grid.setDensitiesFile(file.string()), std::runtime_error);

        // no data
        file = folder / "bad5.csv";

        REQUIRE_THROWS_AS(grid.setDensitiesFile(file.string()), std::runtime_error);

        // non-numeric data
        file = folder / "bad6.csv";

        REQUIRE_THROWS_AS(grid.setDensitiesFile(file.string()), std::runtime_error);
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test translate", "[monad]", Types) {
    TestType grid({2, 3}, {0.5, 1.5});

    // Index=1 → Point=(1,0)
    grid.setDensity(1, 0.5);

    // Point=(0,2)
    grid.translate({1, 2});

    REQUIRE_THAT(grid.getDensity(4), Catch::Matchers::WithinAbs(0.5, NUMERICAL_ZERO));
}

TEST_CASE("monad::Grid2d: Test elementNodes", "[monad]") {
    SECTION("Quad4Grid") {
        const Quad4Grid grid({2, 3}, {0.5, 1.5});

        SECTION("No errors") {
            const Quad4Grid::NodesMatrix expected {
                {0.25, 0.0},
                {0.5, 0.0},
                {0.5, 0.5},
                {0.25, 0.5}
            };

            REQUIRE(grid.elementNodes(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.elementNodes(grid.numElements()), std::out_of_range);
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

            REQUIRE(grid.elementNodes(1).isApprox(expected, NUMERICAL_ZERO));
        }

        SECTION("Invalid index") {
            REQUIRE_THROWS_AS(grid.elementNodes(grid.numElements()), std::out_of_range);
        }
    }
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test measure/area", "[monad]", Types) {
    const TestType grid({2, 3}, {0.5, 1.5});

    REQUIRE_THAT(grid.measure(), Catch::Matchers::WithinAbs(0.75, NUMERICAL_ZERO));
    REQUIRE_THAT(grid.area(), Catch::Matchers::WithinAbs(0.75, NUMERICAL_ZERO));
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test operator==", "[monad]", Types) {
    const TestType grid1({2, 3}, {0.5, 1.5});
    const TestType grid2({2, 3}, {0.5, 1.5});

    REQUIRE(grid1 == grid2);
}

TEMPLATE_LIST_TEST_CASE("monad::Grid2d: Test operator!=", "[monad]", Types) {
    const TestType grid1({2, 3}, {0.5, 1.5});

    SECTION("Different resolution") {
        const TestType grid2({3, 3}, {0.5, 1.5});
        const TestType grid3({2, 4}, {0.5, 1.5});

        REQUIRE(grid1 != grid2);
        REQUIRE(grid1 != grid3);
    }

    SECTION("Different size") {
        const TestType grid2({2, 3}, {0.6, 1.5});
        const TestType grid3({2, 3}, {0.5, 1.6});

        REQUIRE(grid1 != grid2);
        REQUIRE(grid1 != grid3);
    }
}
