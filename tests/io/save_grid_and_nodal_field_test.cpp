#include <tuple>
#include <filesystem>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/io/save_grid_and_nodal_field.hpp"

using namespace monad;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad: Test saveGridAndNodalField", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    
    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    const TestType grid(resolution, size);
    const int numNodes = static_cast<int>(grid.numNodes());

    const std::filesystem::path folder = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path file = folder / "output.msh";

    SECTION("No errors") {
        SECTION("Scalar field") {
            const Eigen::VectorXd field = Eigen::VectorXd::Random(numNodes);

            saveGridAndNodalField(grid, field, file.string());

            REQUIRE(std::filesystem::is_regular_file(file));

            std::filesystem::remove(file);
        }

        SECTION("2D Vector field") {
            const Eigen::MatrixX2d field = Eigen::MatrixX2d::Random(numNodes, 2);

            saveGridAndNodalField(grid, field, file.string());

            REQUIRE(std::filesystem::is_regular_file(file));

            std::filesystem::remove(file);
        }

        SECTION("3D Vector field") {
            const Eigen::MatrixX3d field = Eigen::MatrixX3d::Random(numNodes, 3);

            saveGridAndNodalField(grid, field, file.string());

            REQUIRE(std::filesystem::is_regular_file(file));

            std::filesystem::remove(file);
        }
    }

    SECTION("Invalid field size") {
        const Eigen::VectorXd field = Eigen::VectorXd::Random(numNodes - 1);

        REQUIRE_THROWS_AS(saveGridAndNodalField(grid, field, file.string()), std::invalid_argument);
    }

    SECTION("Invalid file extension") {
        const Eigen::VectorXd field = Eigen::VectorXd::Random(numNodes);

        file = folder / "output.csv";

        REQUIRE_THROWS_AS(saveGridAndNodalField(grid, field, file.string()), std::invalid_argument);
    }

    SECTION("Invalid file path") {
        const Eigen::VectorXd field = Eigen::VectorXd::Random(numNodes);

        file = folder / "NonExistentDir" / "output.msh";

        REQUIRE_THROWS_AS(saveGridAndNodalField(grid, field, file.string()), std::runtime_error);
    }
}
