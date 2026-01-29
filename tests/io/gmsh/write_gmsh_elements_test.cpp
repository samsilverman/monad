#include <sstream>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/io/gmsh/write_gmsh_elements.hpp"

using namespace monad;
using namespace monad::io::gmsh;

TEST_CASE("monad::io::gmsh: Test writeGmshElements", "[monad]") {
    std::ostringstream oss;

    SECTION("Quad4Grid") {
        const Quad4Grid grid({1, 1}, {0.5, 1.0});

        writeGmshElements(oss, grid);

        const std::string expected = "$Elements\n"
                                     "1 1 1 1\n"
                                     "2 1 3 1\n"
                                     "1 1 2 4 3\n"
                                     "$EndElements";

        REQUIRE(oss.str() == expected);
    }

    SECTION("Quad8Grid") {
        const Quad8Grid grid({1, 1}, {0.5, 1.0});

        writeGmshElements(oss, grid);

        const std::string expected = "$Elements\n"
                                     "1 1 1 1\n"
                                     "2 1 16 1\n"
                                     "1 1 2 4 3 5 8 6 7\n"
                                     "$EndElements";

        REQUIRE(oss.str() == expected);
    }

    SECTION("Hex8Grid") {
        const Hex8Grid grid({1, 1, 1}, {0.5, 1.0, 2.0});

        writeGmshElements(oss, grid);

        const std::string expected = "$Elements\n"
                                     "1 1 1 1\n"
                                     "3 1 5 1\n"
                                     "1 1 2 6 5 3 4 8 7\n"
                                     "$EndElements";

        REQUIRE(oss.str() == expected);
    }

    SECTION("Hex20Grid") {
        const Hex20Grid grid({1, 1, 1}, {0.5, 1.0, 2.0});

        writeGmshElements(oss, grid);

        const std::string expected = "$Elements\n"
                                     "1 1 1 1\n"
                                     "3 1 17 1\n"
                                     "1 1 2 6 5 3 4 8 7 9 17 13 18 14 11 16 15 10 19 20 12\n"
                                     "$EndElements";

        REQUIRE(oss.str() == expected);
    }
}
