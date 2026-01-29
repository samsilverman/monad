#include <catch2/catch_test_macros.hpp>
#include "monad/solver/solver_options.hpp"

using namespace monad;

TEST_CASE("monad::SolverOptions: Test defaults", "[monad]") {
    const auto options = SolverOptions::defaults();

    REQUIRE(options.maxIterations == 1000);
    REQUIRE(options.tolerance == 1e-6);
    REQUIRE(options.fields == FieldSave::None);
}

TEST_CASE("monad::SolverOptions: Test operator==", "[monad]") {
    const auto options1 = SolverOptions::defaults();
    const auto options2 = SolverOptions::defaults();

    REQUIRE(options1 == options2);
}

TEST_CASE("monad::SolverOptions: Test operator!=", "[monad]") {
    const auto options1 = SolverOptions::defaults();

    SECTION("Different maxIterations") {
        auto options2 = SolverOptions::defaults();
        options2.maxIterations += 1;

        REQUIRE(options1 != options2);
    }

    SECTION("Different tolerance") {
        auto options2 = SolverOptions::defaults();
        options2.tolerance += 1e-5;

        REQUIRE(options1 != options2);
    }

    SECTION("Different fields") {
        auto options2 = SolverOptions::defaults();
        options2.fields = FieldSave::Total;

        REQUIRE(options1 != options2);
    }
}
