#include <catch2/catch_test_macros.hpp>
#include "monad/grid/hex8_topology.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::grid;

TEST_CASE("monad::grid::Hex8Topology: test numNodes", "[monad]") {
    REQUIRE(Hex8Topology::numNodes({2, 3, 4}) == 60);
}

TEST_CASE("monad::grid::Hex8Topology: test numPeriodicNodes", "[monad]") {
    REQUIRE(Hex8Topology::numPeriodicNodes({2, 3, 4}) == 24);
}

TEST_CASE("monad::grid::Hex8Topology: test node", "[monad]") {
    const Hex8Topology::Point expected(0.25, 0.0, 0.0);
    const Hex8Topology::Point actual = Hex8Topology::node(1, {2, 3, 4}, {0.5, 1.5, 2.0});

    REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::grid::Hex8Topology: test element", "[monad]") {
    const Hex8Topology::ElementList expected{1, 2, 5, 4, 13, 14, 17, 16};
    const Hex8Topology::ElementList actual = Hex8Topology::element(1, {2, 3, 4});

    REQUIRE(actual == expected);
}

TEST_CASE("monad::grid::Hex8Topology: test periodicElement", "[monad]") {
    const Hex8Topology::ElementList expected{1, 0, 2, 3, 7, 6, 8, 9};
    const Hex8Topology::ElementList actual = Hex8Topology::periodicElement(1, {2, 3, 4});

    REQUIRE(actual == expected);
}
