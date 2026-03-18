#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <Eigen/Dense>
#include "monad/detail/constants.hpp"

namespace monad {

    namespace detail {

        /**
         * @brief Symmetrizes a square matrix.
         *
         * Replaces A with ½(A+Aᵀ) to remove numerical asymmetry from
         * floating-point error in computations that should produce a
         * symmetric matrix.
         *
         * @tparam Derived Eigen matrix type.
         *
         * @param[in,out] A Matrix to symmetrize.
         *
         * @throws std::invalid_argument if `A` is not square.
         *
         * @note Use this only to remove numerical noise.
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
         * Uses a Cholesky factorization to test positive
         * definiteness.
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
         * Checks positive semi-definiteness by verifying that all
         * eigenvalues are non-negative up to numerical tolerance.
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
