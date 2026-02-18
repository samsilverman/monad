#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include "monad/detail/constants.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/grid/grid_base.hpp"
#include "monad/solver/solver_options.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"
#include "monad/fem/operator/jacobi_preconditioner.hpp"

namespace monad {

    namespace solver {

        /**
         * @brief Periodic unit cell solver for linear physics problems.
         *
         * This class provides a dimension- and physics-agnostic implementation of a
         * periodic unit cell solver using a matrix-free finite element operator.
         *
         * @tparam Grid Grid class (e.g. Quad4Grid).
         * @tparam Element Element class (e.g. Quad4).
         * @tparam Physics Physics class (e.g. LinearElasticPhysics2d).
         */
        template <class Grid, class Element, class Physics>
        class PeriodicCellSolver {
        public:
            using Material = typename Physics::Material;
            using Kernel = typename Physics::Kernel;
            using OperatorTraits = typename Physics::OperatorTraits;
            using Operator = fem::MatrixFreeOperator<Grid, Element, OperatorTraits>;

            /// @brief Element stiffness matrix type.
            using ElementStiffnessMatrix = typename Kernel::StiffnessMatrix;

            /// @brief Element field matrix type.
            using ElementFieldMatrix = typename Kernel::FieldMatrix;

            /// @brief Global field matrix type.
            using FieldMatrix = typename Physics::FieldMatrix;

            /// @brief Homogenized material tensor type.
            using MaterialTensor = typename Material::MaterialTensor;

            /// @brief Homogenized results type.
            using Results = typename Physics::Results;

            /**
             * @brief Constructs a periodic unit cell solver.
             *
             * @param[in] grid Periodic unit cell grid.
             * @param[in] material Base material model.
             */
            PeriodicCellSolver(const GridBase<Grid, Element> &grid, const Material &material)
                : grid_(grid), material_(material) {
                const auto nodes = grid_.elementNodes(0);

                elementKReference_ = Kernel::lhs(material_, nodes);
                elementFReference_ = Kernel::rhs(material_, nodes);
            }

            /// @brief Periodic unit cell grid_.
            const GridBase<Grid, Element> &grid() const noexcept {
                return grid_;
            }

            /// @brief Base material model.
            const Material &material() const noexcept {
                return material_;
            }

            /**
             * @brief Solves the periodic unit cell problem and computes the homogenized material tensor.
             *
             * The solver computes the total field X decomposed as:
             *
             * X=X̄+X̃
             *
             * - X̄ is the macroscopic field prescribed via uniform macroscopic loading.
             *
             * - X̃ is the microscopic field enforcing equilibrium under periodic boundary conditions.
             *
             * The homogenized material tensor M̄ is then computed via the physics-specific
             * homogenization functional:
             *
             * M̄=1/V∑ₑXₑᵀKₑXₑ
             *
             * @param[in] options Solver options.
             *
             * @returns Homogenized material tensor.
             *
             * @throws std::runtime_error if the linear solver fails to converge.
             */
            Results solve(const SolverOptions &options = SolverOptions::defaults()) const {
                const FieldMatrix XMacro = Physics::macroscopicField(grid_);
                const FieldMatrix XMicro = microscopicField_(options);
                const FieldMatrix X = XMacro + XMicro;

                const MaterialTensor MBar = homogenize_(X);

                return Physics::makeResults(MBar, X, XMacro, XMicro, options.fields);
            }

            /// @brief Equality comparison.
            bool operator==(const PeriodicCellSolver &other) const noexcept {
                return material_ == other.material_ && grid_ == other.grid_;
            }

            /// @brief Inequality comparison.
            bool operator!=(const PeriodicCellSolver &other) const noexcept {
                return !(*this == other);
            }

        private:
            /// @brief Periodic unit cell grid_.
            GridBase<Grid, Element> grid_;

            /// @brief Base material model.
            Material material_;

            /// @brief Reference element stiffness matrix for unit density.
            ElementStiffnessMatrix elementKReference_;

            /// @brief Reference element source matrix for unit density.
            ElementFieldMatrix elementFReference_;

            /**
             * @brief Builds the reduced right-hand side source matrix.
             *
             * Constructs the source matrix F in the reduced linear system
             *
             * KX̃=F
             *
             * where K acts only on a reduced set of unconstrained periodic dofs.
             *
             * @param[in] K Matrix-free stiffness operator.
             *
             * @returns Reduced right-hand side source matrix.
             */
            FieldMatrix buildReducedRhs_(const Operator &K) const noexcept {
                const std::size_t numPeriodicNodes = grid_.numPeriodicNodes();
                const std::size_t numPeriodicDofs = OperatorTraits::NumNodeDofs * numPeriodicNodes;
                const std::size_t numReducedDofs = numPeriodicDofs - OperatorTraits::NumFixedDofs;

                FieldMatrix reducedF = FieldMatrix::Zero(static_cast<int>(numReducedDofs), static_cast<int>(Physics::NumMacroFields));

                const auto &elementDofs = K.elementDofs();

                for (std::size_t i = 0; i < grid_.numElements(); ++i) {
                    const double density = grid_.getDensity(i);
                    const ElementFieldMatrix elementF = density * elementFReference_;

                    const auto &dofs = elementDofs[i];

                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const int dof = dofs[j];

                        if (dof >= 0) {
                            reducedF(dof, Eigen::indexing::all) += elementF.row(static_cast<int>(j));
                        }
                    }
                }

                return reducedF;
            }

            /**
             * @brief Builds the expanded left-hand side microscopic field.
             *
             * Expands the microscopic field matrix X̃ in the reduced linear system
             *
             * KX̃=F
             *
             * from a reduced set of unconstrained periodic dofs to global dofs.
             *
             * @param[in] reducedX Reduced microscopic field.
             *
             * @returns Expanded left-hand side microscopic field.
             */
            FieldMatrix buildExpandedLhs_(const FieldMatrix &reducedX) const noexcept {
                // Expand reduced to periodic
                const std::size_t numPeriodicNodes = grid_.numPeriodicNodes();
                const std::size_t numPeriodicDofs = OperatorTraits::NumNodeDofs * numPeriodicNodes;
                const std::size_t numReducedDofs = numPeriodicDofs - OperatorTraits::NumFixedDofs;

                FieldMatrix periodicX = FieldMatrix::Zero(static_cast<int>(numPeriodicDofs), static_cast<int>(Physics::NumMacroFields));

                for (std::size_t i = 0; i < numReducedDofs; ++i) {
                    const std::size_t iExpanded = OperatorTraits::expandedDof(i, numPeriodicNodes);
                    periodicX.row(static_cast<int>(iExpanded)) = reducedX.row(static_cast<int>(i));
                }

                // Expand periodic to global
                const std::size_t numNodes = grid_.numNodes();
                const std::size_t numDofs = OperatorTraits::NumNodeDofs * numNodes;

                FieldMatrix X = FieldMatrix::Zero(static_cast<int>(numDofs), static_cast<int>(Physics::NumMacroFields));

                for (std::size_t i = 0; i < grid_.numElements(); ++i) {
                    const auto element = grid_.element(i);
                    const auto periodicElement = grid_.periodicElement(i);

                    const auto dofs = OperatorTraits::dofs(element, numNodes);
                    const auto periodicDofs = OperatorTraits::dofs(periodicElement, numPeriodicNodes);

                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const std::size_t dof = dofs[j];
                        const std::size_t periodicDof = periodicDofs[j];

                        X.row(static_cast<int>(dof)) = periodicX.row(static_cast<int>(periodicDof));
                    }
                }

                return X;
            }

            /**
             * @brief Microscopic correction field.
             *
             * Solves the reduced linear system:
             *
             * KX̃=F
             *
             * using a matrix-free solver then expands
             * the solution to global dofs.
             *
             * @param[in] options Solver options.
             *
             * @returns Microscopic correction field.
             *
             * @throws std::runtime_error if the solver fails to converge.
             *
             * @note If the global stiffness matric is symmetric, conjugate gradient
             * is used, otherwise biconjugate gradient stabilized (BiCSSTAB) is used.
             */
            FieldMatrix microscopicField_(const SolverOptions &options) const {
                const Operator K(grid_, elementKReference_);
                const FieldMatrix reducedF = buildReducedRhs_(K);

                using preconditioner = fem::JacobiPreconditioner<Operator>;

                Eigen::ConjugateGradient<Operator, Eigen::Lower | Eigen::Upper, preconditioner> cg;
                cg.setMaxIterations(options.maxIterations);
                cg.setTolerance(options.tolerance);
                cg.compute(K);

                FieldMatrix reducedX = cg.solve(reducedF);

                if (cg.info() != Eigen::Success) {
                    std::string error_message = "Solver failed to converge.";
                    if (cg.info() == Eigen::NoConvergence) {
                        error_message += " No Convergence - Max iterations reached or tolerance not met.";
                    }
                    else if (cg.info() == Eigen::NumericalIssue) {
                        error_message += " Numerical Issue - Solver encountered stability problems.";
                    }

                    throw std::runtime_error(error_message);
                }

                return buildExpandedLhs_(reducedX);
            }

            /**
             * @brief Homogenized material tensor.
             *
             * The homogenized material tensor M̄ is obtained via the Hill-Mandel lemma,
             *
             * M̄=1/V∑ₑXₑᵀKₑXₑ
             *
             * where the total field X is the sum of the macroscopic and microscopic components, X=X̄+X̃.
             *
             * @param[in] X Nodal total field.
             *
             * @returns Homogenized material tensor.
             */
            MaterialTensor homogenize_(const FieldMatrix &X) const {
                MaterialTensor MBar = MaterialTensor::Zero();

                const std::size_t numNodes = grid_.numNodes();

                for (std::size_t i = 0; i < grid_.numElements(); ++i) {
                    const double density = grid_.getDensity(i);
                    const auto elementK = density * elementKReference_;
                    const auto element = grid_.element(i);
                    const auto dofs = OperatorTraits::dofs(element, numNodes);
                    const auto elementX = X(detail::arrayToEigen(dofs), Eigen::indexing::all);

                    MBar.noalias() += elementX.transpose() * elementK * elementX;
                }

                MBar /= grid_.measure();

                // Remove numerical artifacts
                detail::symmetrize(MBar);

                return MBar;
            }
        };

    } // namespace solver

} // namespace monad
