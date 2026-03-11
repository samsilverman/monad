#include <catch2/catch_test_macros.hpp>
#include "monad/grid/hex20_topology.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::grid;

TEST_CASE("monad::grid::Hex20Topology: test numNodes", "[monad]") {
    REQUIRE(Hex20Topology::numNodes({2, 3, 4}) == 193);
}

TEST_CASE("monad::grid::Hex20Topology: test numPeriodicNodes", "[monad]") {
    REQUIRE(Hex20Topology::numPeriodicNodes({2, 3, 4}) == 96);
}

TEST_CASE("monad::grid::Hex20Topology: test node", "[monad]") {
    const Hex20Topology::Point expected(0.25, 0.0, 0.0);
    const Hex20Topology::Point actual = Hex20Topology::node(1, {2, 3, 4}, {0.5, 1.5, 2.0});

    REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::grid::Hex20Topology: test element", "[monad]") {
    const Hex20Topology::ElementList expected{1, 2, 5, 4, 13, 14, 17, 16, 61, 102, 63, 101, 69, 111, 71, 110, 146, 147, 150, 149};
    const Hex20Topology::ElementList actual = Hex20Topology::element(1, {2, 3, 4});

    REQUIRE(actual == expected);
}

TEST_CASE("monad::grid::Hex20Topology: test periodicElement", "[monad]") {
    const Hex20Topology::ElementList expected{1, 0, 2, 3, 7, 6, 8, 9, 25, 48, 27, 49, 31, 54, 33, 55, 73, 72, 74, 75};
    const Hex20Topology::ElementList actual = Hex20Topology::periodicElement(1, {2, 3, 4});

    REQUIRE(actual == expected);
}
