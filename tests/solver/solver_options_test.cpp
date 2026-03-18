#include <catch2/catch_test_macros.hpp>
#include "monad/solver/solver_options.hpp"

using namespace monad;

TEST_CASE("monad::SolverOptions: Test defaults", "[monad]") {
    const SolverOptions options;

    REQUIRE(options.maxIterations == 1000);
    REQUIRE(options.tolerance == 1e-6);
    REQUIRE(options.fields == FieldSave::None);
}

TEST_CASE("monad::SolverOptions: Test operator==", "[monad]") {
    const SolverOptions options1;
    const SolverOptions options2;

    REQUIRE(options1 == options2);
}

TEST_CASE("monad::SolverOptions: Test operator!=", "[monad]") {
    const SolverOptions options1;
    SolverOptions options2;

    SECTION("Different maxIterations") {
        options2.maxIterations += 1;

        REQUIRE(options1 != options2);
    }

    SECTION("Different tolerance") {
        options2.tolerance += 1e-5;

        REQUIRE(options1 != options2);
    }

    SECTION("Different fields") {
        options2.fields = FieldSave::Total;

        REQUIRE(options1 != options2);
    }
}
