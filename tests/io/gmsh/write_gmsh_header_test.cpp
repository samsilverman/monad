#include <sstream>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "monad/io/gmsh/write_gmsh_header.hpp"

using namespace monad::io::gmsh;

TEST_CASE("monad::io::gmsh: Test writeGmshHeader", "[monad]") {
    std::ostringstream oss;

    writeGmshHeader(oss);

    const std::string expected = "$MeshFormat\n"
                                 "4.1 0 8\n"
                                 "$EndMeshFormat";

    REQUIRE(oss.str() == expected);
}
