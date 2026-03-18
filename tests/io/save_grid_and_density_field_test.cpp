#include <tuple>
#include <filesystem>
#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/field/density_field.hpp"
#include "monad/io/save_grid_and_density_field.hpp"

using namespace monad;
using namespace monad::field;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad: Test saveGridAndDensityField", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    
    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    const TestType grid(resolution, size);
    const DensityField<TestType::Dim> densityField(resolution);

    const std::filesystem::path folder = std::filesystem::path(__FILE__).parent_path();
    std::filesystem::path file = folder / "output.msh";

    SECTION("No errors") {
        saveGridAndDensityField(grid, densityField, file.string());

        REQUIRE(std::filesystem::is_regular_file(file));

        std::filesystem::remove(file);
    }

    SECTION("Invalid field size") {
        resolution[0] += 1;

        REQUIRE_THROWS_AS(saveGridAndDensityField(grid, DensityField<TestType::Dim>(resolution), file.string()), std::invalid_argument);
    }

    SECTION("Invalid file extension") {
        file = folder / "output.csv";

        REQUIRE_THROWS_AS(saveGridAndDensityField(grid, densityField, file.string()), std::invalid_argument);
    }

    SECTION("Invalid file path") {
        file = folder / "NonExistentDir" / "output.msh";

        REQUIRE_THROWS_AS(saveGridAndDensityField(grid, densityField, file.string()), std::runtime_error);
    }
}
