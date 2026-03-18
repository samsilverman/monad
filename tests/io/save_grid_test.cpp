#include <tuple>
#include <filesystem>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/io/save_grid.hpp"

using namespace monad;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad: Test saveGrid", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    
    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);

    const std::filesystem::path folder = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path file = folder / "output.msh";

    SECTION("No errors") {
        saveGrid(grid, file.string());

        REQUIRE(std::filesystem::is_regular_file(file));

        std::filesystem::remove(file);
    }

    SECTION("Invalid file extension") {
        file = folder / "output.csv";

        REQUIRE_THROWS_AS(saveGrid(grid, file.string()), std::invalid_argument);
    }

    SECTION("Invalid file path") {
        file = folder / "NonExistentDir" / "output.msh";

        REQUIRE_THROWS_AS(saveGrid(grid, file.string()), std::runtime_error);
    }
}
