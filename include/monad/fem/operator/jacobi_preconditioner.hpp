#pragma once

#include <Eigen/Core>

namespace monad {

    namespace fem {

        /**
         * @brief Jacobi (diagonal) preconditioner for the matrix-free operator.
         *
         * This preconditioner approximates the inverse of the global stiffness matrix A
         * by inverting only its diagonal:
         *
         * A⁻¹≈diag(A)⁻¹
         *
         * @tparam Operator Matrix-free operator class (e.g. MatrixFreeOperator<Quad4Grid, Quad4, LinearElasticMatrixFreeOperatorTraits<2>>).
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
             * @brief Computes the diagonal of the global stiffness matrix.
             *
             * @param[in] op Matrix-free operator.
             *
             * @returns Reference to this preconditioner.
             */
            JacobiPreconditioner &compute(const Operator &op) noexcept {
                const int n = op.rows();
                diagonal_.setZero(n);

                const auto &elementKReference = op.elementKReference();
                const auto &elementDofs = op.elementDofs();
                const auto &densities = op.densities();

                // Loop over all elements and accumulate their diagonal contributions
                for (std::size_t i = 0; i < elementDofs.size(); ++i) {
                    const double density = densities[i];
                    const auto &dofs = elementDofs[i];

                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const int globalDof = dofs[j];
                        const int localDof = static_cast<int>(j);

                        if (globalDof >= 0) {
                            diagonal_(globalDof) += density * elementKReference(localDof, localDof);
                        }
                    }
                }

                return *this;
            }

            /**
             * @brief Applies the Jacobi preconditioner.
             *
             * Computes
             *
             * x=A⁻¹b≈diag(A)⁻¹b
             *
             * via element-wise division.
             *
             * @tparam Rhs Type of the right-hand side vector.
             *
             * @param[in] b Right-hand side vector.
             *
             * @returns Preconditioned vector.
             */
            template <typename Rhs>
            inline Vector solve(const Eigen::MatrixBase<Rhs> &b) const noexcept {
                return b.cwiseQuotient(diagonal_);
            }

            /// @brief Status of the preconditioner.
            Eigen::ComputationInfo info() const noexcept {
                return Eigen::Success;
            }

        private:
            /// @brief Diagonal of the global stiffness matrix.
            Vector diagonal_;
        };

    } // namespace fem

} // namespace monad
