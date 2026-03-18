#include <tuple>
#include <stdexcept>
#include <filesystem>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "monad/field/field_aliases.hpp"
#include "monad/field/make_density_field_from_csv.hpp"

using namespace monad;

TEST_CASE("monad: Test makeDensityFieldFromCsv", "[monad]") {
    using DensityList = typename DensityField2d::DensityList;

    DensityField2d densityField({2, 3});

    const auto folder = std::filesystem::path(__FILE__).parent_path() / "csv";

    SECTION("No errors") {
        const auto file = folder / "good.csv";

        auto densityField = makeDensityFieldFromCsv(file.string());

        const DensityList expected {
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6
        };

        REQUIRE_THAT(densityField.densities(), Catch::Matchers::Approx(expected).margin(NUMERICAL_ZERO));
    }

    SECTION("Invalid csv data") {
        // data < 0
        auto file = folder / "bad1.csv";

        REQUIRE_THROWS_AS(makeDensityFieldFromCsv(file.string()), std::runtime_error);

        // data > 1
        file = folder / "bad2.csv";

        REQUIRE_THROWS_AS(makeDensityFieldFromCsv(file.string()), std::runtime_error);

        // not enough data
        file = folder / "bad3.csv";

        REQUIRE_THROWS_AS(makeDensityFieldFromCsv(file.string()), std::runtime_error);

        // too much data
        file = folder / "bad4.csv";

        REQUIRE_THROWS_AS(makeDensityFieldFromCsv(file.string()), std::runtime_error);

        // no data
        file = folder / "bad5.csv";

        REQUIRE_THROWS_AS(makeDensityFieldFromCsv(file.string()), std::runtime_error);

        // non-numeric data
        file = folder / "bad6.csv";

        REQUIRE_THROWS_AS(makeDensityFieldFromCsv(file.string()), std::runtime_error);
    }
}
