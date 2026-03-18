#pragma once

#include <Eigen/Core>

namespace monad {

    namespace solver {

        /**
         * @brief Jacobi preconditioner for a matrix-free stiffness operator.
         *
         * This preconditioner approximates the inverse of the reduced global
         * stiffness operator A by inverting only its diagonal:
         *
         * ```text
         * A⁻¹≈diag(A)⁻¹
         * ```
         *
         * @tparam Operator Matrix-free operator type (e.g. MatrixFreeOperator<Quad4Grid, LinearElasticDofTraits<2>>).
         */
        template <class Operator>
        class JacobiPreconditioner {
        public:
            using Scalar = typename Operator::Scalar;
            using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

            /// @brief Default constructor.
            JacobiPreconditioner() = default;

            /**
             * @brief Constructs the Jacobi preconditioner.
             *
             * @param[in] op Matrix-free operator.
             */
            explicit JacobiPreconditioner(const Operator &op) noexcept {
                compute(op);
            }

            /**
             * @brief Computes the diagonal preconditioner from a matrix-free operator.
             *
             * @param[in] op Matrix-free operator.
             *
             * @returns Reference to this preconditioner.
             */
            JacobiPreconditioner &compute(const Operator &op) noexcept {
                const int n = op.rows();
                diagonal_.setZero(n);

                const auto &elementKReference = op.elementKReference();
                const auto &dofMap = op.dofMap();
                const auto &densityField = op.densityField();

                // Assemble the reduced global diagonal by accumulating the
                // density-scaled diagonal entries of each element matrix
                for (std::size_t i = 0; i < dofMap.size(); ++i) {
                    const double density = densityField.getDensity(i);
                    const auto &reducedDofs = dofMap.reducedDofs(i);

                    for (std::size_t j = 0; j < reducedDofs.size(); ++j) {
                        const int reducedDof = reducedDofs[j];
                        const int localReducedDof = static_cast<int>(j);

                        if (reducedDof >= 0) {
                            diagonal_(reducedDof) += density * elementKReference(localReducedDof, localReducedDof);
                        }
                    }
                }

                // Jacobi requires invertible diagonal entries
                for (int i = 0; i < diagonal_.size(); ++i) {
                    if (diagonal_(i) == Scalar(0)) {
                        info_ = Eigen::NumericalIssue;
                        break;
                    }
                }

                return *this;
            }

            /**
             * @brief Applies the Jacobi preconditioner.
             *
             * Computes
             *
             * ```text
             * x=A⁻¹b≈diag(A)⁻¹b
             * ```
             *
             * by element-wise division.
             *
             * @tparam Rhs Type of the right-hand side vector.
             *
             * @param[in] b Right-hand side vector.
             *
             * @returns Preconditioned vector.
             */
            template <typename Rhs>
            inline Vector solve(const Eigen::MatrixBase<Rhs> &b) const {
                return b.cwiseQuotient(diagonal_);
            }

            /// @brief Status of the preconditioner.
            Eigen::ComputationInfo info() const noexcept {
                return info_;
            }

        private:
            /// @brief Diagonal of the global stiffness matrix.
            Vector diagonal_;

            /// @brief Status of the preconditioner.
            Eigen::ComputationInfo info_ = Eigen::Success;
        };

    } // namespace solver

} // namespace monad
