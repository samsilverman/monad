#include <sstream>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/io/gmsh/write_gmsh_nodes.hpp"

using namespace monad;
using namespace monad::io::gmsh;

TEST_CASE("monad::io::gmsh: Test writeGmshNodes", "[monad]") {
    std::ostringstream oss;

    SECTION("Quad4Grid") {
        const Quad4Grid grid({1, 1}, {0.5, 1.0});

        writeGmshNodes(oss, grid);

        const std::string expected = "$Nodes\n"
                                     "1 4 1 4\n"
                                     "2 1 0 4\n"
                                     "1\n"
                                     "2\n"
                                     "3\n"
                                     "4\n"
                                     "0 0 0\n"
                                     "0.5 0 0\n"
                                     "0 1 0\n"
                                     "0.5 1 0\n"
                                     "$EndNodes";

        REQUIRE(oss.str() == expected);
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({1, 1}, {0.5, 1.0});

        writeGmshNodes(oss, grid);

        const std::string expected = "$Nodes\n"
                                     "1 8 1 8\n"
                                     "2 1 0 8\n"
                                     "1\n"
                                     "2\n"
                                     "3\n"
                                     "4\n"
                                     "5\n"
                                     "6\n"
                                     "7\n"
                                     "8\n"
                                     "0 0 0\n"
                                     "0.5 0 0\n"
                                     "0 1 0\n"
                                     "0.5 1 0\n"
                                     "0.25 0 0\n"
                                     "0.25 1 0\n"
                                     "0 0.5 0\n"
                                     "0.5 0.5 0\n"
                                     "$EndNodes";

        REQUIRE(oss.str() == expected);
    }

    SECTION("Hex8Grid") {
        const Hex8Grid grid({1, 1, 1}, {0.5, 1.0, 2.0});

        writeGmshNodes(oss, grid);

        const std::string expected = "$Nodes\n"
                                     "1 8 1 8\n"
                                     "3 1 0 8\n"
                                     "1\n"
                                     "2\n"
                                     "3\n"
                                     "4\n"
                                     "5\n"
                                     "6\n"
                                     "7\n"
                                     "8\n"
                                     "0 0 0\n"
                                     "0.5 0 0\n"
                                     "0 1 0\n"
                                     "0.5 1 0\n"
                                     "0 0 2\n"
                                     "0.5 0 2\n"
                                     "0 1 2\n"
                                     "0.5 1 2\n"
                                     "$EndNodes";

        REQUIRE(oss.str() == expected);
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({1, 1, 1}, {0.5, 1.0, 2.0});

        writeGmshNodes(oss, grid);

        const std::string expected = "$Nodes\n"
                                     "1 20 1 20\n"
                                     "3 1 0 20\n"
                                     "1\n"
                                     "2\n"
                                     "3\n"
                                     "4\n"
                                     "5\n"
                                     "6\n"
                                     "7\n"
                                     "8\n"
                                     "9\n"
                                     "10\n"
                                     "11\n"
                                     "12\n"
                                     "13\n"
                                     "14\n"
                                     "15\n"
                                     "16\n"
                                     "17\n"
                                     "18\n"
                                     "19\n"
                                     "20\n"
                                     "0 0 0\n"
                                     "0.5 0 0\n"
                                     "0 1 0\n"
                                     "0.5 1 0\n"
                                     "0 0 2\n"
                                     "0.5 0 2\n"
                                     "0 1 2\n"
                                     "0.5 1 2\n"
                                     "0.25 0 0\n"
                                     "0.25 1 0\n"
                                     "0.25 0 2\n"
                                     "0.25 1 2\n"
                                     "0 0.5 0\n"
                                     "0.5 0.5 0\n"
                                     "0 0.5 2\n"
                                     "0.5 0.5 2\n"
                                     "0 0 1\n"
                                     "0.5 0 1\n"
                                     "0 1 1\n"
                                     "0.5 1 1\n"
                                     "$EndNodes";

        REQUIRE(oss.str() == expected);
    }
}
