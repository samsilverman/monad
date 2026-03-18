#include <tuple>
#include <sstream>
#include <string>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/field/field_aliases.hpp"
#include "monad/io/gmsh/write_gmsh_densities.hpp"

using namespace monad;
using namespace monad::io::gmsh;

using Types = std::tuple<DensityField2d, DensityField3d>;

TEMPLATE_LIST_TEST_CASE("monad::io::gmsh: Test writeGmshDensities", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;

    Resolution resolution;
    resolution.fill(1);

    TestType densityField(resolution);
    densityField.setRandom(1234);

    std::ostringstream oss;

    writeGmshDensities(oss, densityField);

    const std::string expected = "$ElementData\n"
                                 "1\n"
                                 "\"Density\"\n"
                                 "0\n"
                                 "3\n"
                                 "0\n"
                                 "1\n"
                                 "1\n"
                                 "1 0.191519\n"
                                 "$EndElementData";

    REQUIRE(oss.str() == expected);
}
