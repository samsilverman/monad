#include <array>
#include <cstddef>
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include "monad/detail/eigen_utils.hpp"

using namespace monad::detail;

TEST_CASE("monad::detail: Test arrayToEigen", "[monad]") {
    const std::array<std::size_t, 3> array{1, 2, 3};

    const Eigen::Vector3i expected(1, 2, 3);

    REQUIRE(arrayToEigen(array) == expected);
}

TEST_CASE("monad::detail: Test symmetrize/isSymmetric", "[monad]") {
    Eigen::Matrix3d A = Eigen::Matrix3d::Random();

    REQUIRE(!isSymmetric(A));

    symmetrize(A);

    REQUIRE(isSymmetric(A));
}

TEST_CASE("monad::detail: Test isPD/isPSD", "[monad]") {
    Eigen::Matrix3d A;
    A << 1, 2, 3,
         4, 5, 6,
         0, 0, 0;

    REQUIRE(!isPD(A));
    REQUIRE(!isPSD(A));

    // Make PSD
    A = A.transpose() * A;

    REQUIRE(!isPD(A));
    REQUIRE(isPSD(A));

    // Make PD
    A += Eigen::Matrix3d::Identity();

    REQUIRE(isPD(A));
    REQUIRE(isPSD(A));
}
