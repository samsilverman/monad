#include <tuple>
#include <sstream>
#include <string>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/io/gmsh/write_gmsh_densities.hpp"

using namespace monad;
using namespace monad::io::gmsh;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::io::gmsh: Test writeGmshDensities", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;

    Resolution resolution;
    resolution.fill(1);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);
    grid.setDensitiesRandom(1234);

    std::ostringstream oss;

    writeGmshDensities(oss, grid);

    const std::string expected = "$ElementData\n"
                                 "1\n"
                                 "\"Density\"\n"
                                 "0\n"
                                 "3\n"
                                 "0\n"
                                 "1\n"
                                 "1\n"
                                 "1 0.497664\n"
                                 "$EndElementData";

    REQUIRE(oss.str() == expected);
}
