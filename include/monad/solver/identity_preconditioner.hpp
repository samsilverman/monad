#pragma once

#include <Eigen/Sparse>

namespace monad {

    namespace solver {

        /**
         * @brief Identity preconditioner for a matrix-free stiffness operator.
         *
         * This preconditioner approximates the inverse of the reduced global
         * stiffness operator A by the identity matrix:
         *
         * ```text
         * A⁻¹≈I
         * ```
         */
        using IdentityPreconditioner = Eigen::IdentityPreconditioner;

    } // namespace solver

} // namespace monad
