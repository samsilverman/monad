#include <catch2/catch_test_macros.hpp>
#include "monad/grid/quad4_topology.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::grid;

TEST_CASE("monad::grid::Quad4Topology: test numNodes", "[monad]") {
    REQUIRE(Quad4Topology::numNodes({2, 3}) == 12);
}

TEST_CASE("monad::grid::Quad4Topology: test numPeriodicNodes", "[monad]") {
    REQUIRE(Quad4Topology::numPeriodicNodes({2, 3}) == 6);
}

TEST_CASE("monad::grid::Quad4Topology: test node", "[monad]") {
    const Quad4Topology::Point expected(0.25, 0.0);
    const Quad4Topology::Point actual = Quad4Topology::node(1, {2, 3}, {0.5, 1.5});

    REQUIRE(actual.isApprox(expected, NUMERICAL_ZERO));
}

TEST_CASE("monad::grid::Quad4Topology: test element", "[monad]") {
    const Quad4Topology::ElementList expected{1, 2, 5, 4};
    const Quad4Topology::ElementList actual = Quad4Topology::element(1, {2, 3});

    REQUIRE(actual == expected);
}

TEST_CASE("monad::grid::Quad4Topology: test periodicElement", "[monad]") {
    const Quad4Topology::ElementList expected{1, 0, 2, 3};
    const Quad4Topology::ElementList actual = Quad4Topology::periodicElement(1, {2, 3});

    REQUIRE(actual == expected);
}
