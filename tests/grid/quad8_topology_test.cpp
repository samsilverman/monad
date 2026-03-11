#include <catch2/catch_test_macros.hpp>
#include "monad/grid/quad8_topology.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::grid;

TEST_CASE("monad::grid::Quad8Topology: test numNodes", "[monad]") {
    REQUIRE(Quad8Topology::numNodes({2, 3}) == 29);
}

TEST_CASE("monad::grid::Quad8Topology: test numPeriodicNodes", "[monad]") {
    REQUIRE(Quad8Topology::numPeriodicNodes({2, 3}) == 18);
}

TEST_CASE("monad::grid::Quad8Topology: test node", "[monad]") {
    const Quad8Topology::Point expected(0.25, 0.0);
    const Quad8Topology::Point actual = Quad8Topology::node(1, {2, 3}, {0.5, 1.5});

    REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::grid::Quad8Topology: test element", "[monad]") {
    const Quad8Topology::ElementList expected{1, 2, 5, 4, 13, 22, 15, 21};
    const Quad8Topology::ElementList actual = Quad8Topology::element(1, {2, 3});

    REQUIRE(actual == expected);
}

TEST_CASE("monad::grid::Quad8Topology: test periodicElement", "[monad]") {
    const Quad8Topology::ElementList expected{1, 0, 2, 3, 7, 12, 9, 13};
    const Quad8Topology::ElementList actual = Quad8Topology::periodicElement(1, {2, 3});

    REQUIRE(actual == expected);
}
