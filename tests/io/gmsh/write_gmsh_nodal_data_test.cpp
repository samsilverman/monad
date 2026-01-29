#include <sstream>
#include <string>
#include <stdexcept>
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include "monad/io/gmsh/write_gmsh_nodal_field.hpp"

using namespace monad::io::gmsh;

TEST_CASE("monad::io::gmsh: Test writeGmshNodalFields", "[monad]") {
    std::ostringstream oss;

    SECTION("Scalar field") {
        const Eigen::Vector3d field(0.1, 0.2, 0.3);

        writeGmshNodalField(oss, field, "name");

        const std::string expected = "$NodeData\n"
                                     "1\n"
                                     "\"name\"\n"
                                     "0\n"
                                     "3\n"
                                     "0\n"
                                     "1\n"
                                     "3\n"
                                     "1 0.1\n"
                                     "2 0.2\n"
                                     "3 0.3\n"
                                     "$EndNodeData";

        REQUIRE(oss.str() == expected);
    }

    SECTION("2D vector field") {
        const Eigen::Matrix2d field {
            {0.1, 0.2},
            {0.3, 0.4},
        };

        writeGmshNodalField(oss, field);

        const std::string expected = "$NodeData\n"
                                     "0\n"
                                     "0\n"
                                     "3\n"
                                     "0\n"
                                     "3\n"
                                     "2\n"
                                     "1 0.1 0.2 0\n"
                                     "2 0.3 0.4 0\n"
                                     "$EndNodeData";

        REQUIRE(oss.str() == expected);
    }

    SECTION("3D vector field") {
        const Eigen::Matrix3d field {
            {0.1, 0.2, 0.3},
            {0.4, 0.5, 0.6},
            {0.7, 0.8, 0.9}
        };

        writeGmshNodalField(oss, field, "name");

        const std::string expected = "$NodeData\n"
                                     "1\n"
                                     "\"name\"\n"
                                     "0\n"
                                     "3\n"
                                     "0\n"
                                     "3\n"
                                     "3\n"
                                     "1 0.1 0.2 0.3\n"
                                     "2 0.4 0.5 0.6\n"
                                     "3 0.7 0.8 0.9\n"
                                     "$EndNodeData";

        REQUIRE(oss.str() == expected);
    }

    SECTION("invalid vector field") {
        const Eigen::Matrix4d field = Eigen::Matrix4d::Random();

        REQUIRE_THROWS_AS(writeGmshNodalField(oss, field, "name"), std::invalid_argument);
    }
}
