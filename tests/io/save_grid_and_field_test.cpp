#include <tuple>
#include <filesystem>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/io/save_grid_and_field.hpp"

using namespace monad;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad: Test saveGridAndField", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    
    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);
    grid.setDensitiesRandom();
    const int numNodes = static_cast<int>(grid.numNodes());

    const std::filesystem::path folder = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path file = folder / "output.msh";

    SECTION("No errors") {
        const Eigen::MatrixX2d field = Eigen::MatrixX2d::Random(numNodes, 2);

        saveGridAndField(grid, field, file.string());

        REQUIRE(std::filesystem::is_regular_file(file));

        std::filesystem::remove(file);
    }

    SECTION("Invalid field size") {
        const Eigen::MatrixX4d field = Eigen::MatrixX4d::Random(numNodes, 4);

        REQUIRE_THROWS_AS(saveGridAndField(grid, field, file.string()), std::invalid_argument);
    }

    SECTION("Invalid file extension") {
        const Eigen::MatrixX2d field = Eigen::MatrixX2d::Random(numNodes, 2);

        file = folder / "output.csv";

        REQUIRE_THROWS_AS(saveGridAndField(grid, field, file.string()), std::invalid_argument);
    }

    SECTION("Invalid file path") {
        const Eigen::MatrixX2d field = Eigen::MatrixX2d::Random(numNodes, 2);

        file = folder / "NonExistentDir" / "output.msh";

        REQUIRE_THROWS_AS(saveGridAndField(grid, field, file.string()), std::runtime_error);
    }
}
