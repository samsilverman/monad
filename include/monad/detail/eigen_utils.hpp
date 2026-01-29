#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <Eigen/Dense>
#include "monad/detail/constants.hpp"

namespace monad {

    namespace detail {

        /**
         * @brief Converts an array to an Eigen vector.
         *
         * @tparam N Size of array.
         *
         * @param[in] array Array.
         *
         * @returns Eigen vector of `array`.
         */
        template <std::size_t N>
        inline Eigen::Vector<int, N> arrayToEigen(const std::array<std::size_t, N> &array) noexcept {
            Eigen::Vector<int, N> vector;

            for (std::size_t i = 0; i < N; ++i) {
                vector[static_cast<int>(i)] = static_cast<int>(array[i]);
            }

            return vector;
        }

        /**
         * @brief Symmetrizes a square matrix.
         *
         * This function replaces A with ½(A+Aᵀ), removing numerical
         * asymmetry introduced by floating-point error in computations
         * that theoretically produce symmetric matrices.
         *
         * @tparam Derived Eigen matrix type.
         *
         * @param[in,out] A Matrix to symmetrize.
         *
         * @throws std::invalid_argument if `A` is not square.
         *
         * @note Use this only to clean up numerical noise.
         */
        template <typename Derived>
        inline void symmetrize(Eigen::MatrixBase<Derived> &A) {
            if (A.rows() != A.cols()) {
                throw std::invalid_argument("A must be a square matrix.");
            }

            A = 0.5 * (A + A.transpose()).eval();
        }

        /**
         * @brief Checks if a matrix is symmetric.
         *
         * @tparam Derived Eigen matrix type.
         *
         * @param[in] A Matrix to check.
         *
         * @returns `true` if `A` is symmetric, `false` otherwise.
         */
        template <typename Derived>
        inline bool isSymmetric(const Eigen::MatrixBase<Derived> &A) noexcept {
            if (A.rows() != A.cols()) {
                return false;
            }

            return A.isApprox(A.transpose(), NUMERICAL_ZERO);
        }

        /**
         * @brief Checks if a matrix is positive definite.
         *
         * This function checks if a matrix is positive definite
         * by attempting to perform a Cholesky factorization.
         *
         * @tparam Derived Eigen matrix type.
         *
         * @param[in] A Matrix to check.
         *
         * @returns `true` if `A` is is positive definite, `false` otherwise.
         */
        template <typename Derived>
        inline bool isPD(const Eigen::MatrixBase<Derived> &A) noexcept {
            if (!isSymmetric(A)) {
                return false;
            }

            using PlainObject = typename Derived::PlainObject;
            Eigen::LLT<PlainObject> llt(A.eval());

            return llt.info() == Eigen::Success;
        }

        /**
         * @brief Checks if a matrix is positive semi-definite.
         *
         * This function checks if a matrix is positive semi-definite
         * by verifying that all its eigenvalues are non-negative.
         *
         * @tparam Derived Eigen matrix type.
         *
         * @param[in] A Matrix to check.
         *
         * @returns `true` if `A` is is positive semi-definite, `false` otherwise.
         */
        template <typename Derived>
        inline bool isPSD(const Eigen::MatrixBase<Derived> &A) noexcept {
            if (!isSymmetric(A)) {
                return false;
            }

            using PlainObject = typename Derived::PlainObject;
            Eigen::SelfAdjointEigenSolver<PlainObject> solver(A.eval());

            if (solver.info() != Eigen::Success) {
                return false;
            }

            return (solver.eigenvalues().array() >= -NUMERICAL_ZERO).all();
        }

    } // namespace detail

} // namespace monad
