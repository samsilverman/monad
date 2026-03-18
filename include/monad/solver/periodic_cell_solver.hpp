#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include "monad/field/density_field.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/solver/solver_options.hpp"
#include "monad/solver/jacobi_preconditioner.hpp"

namespace monad {

    namespace solver {

        /**
         * @brief Stateless periodic-cell solver for structured linear homogenization problems.
         *
         * @tparam Grid Grid type (e.g. Quad4Grid).
         * @tparam Policy Physics policy type (e.g. LinearElasticPolicy<Quad4>).
         */
        template <class Grid, class Policy>
        class PeriodicCellSolver {
        public:
            using Material = typename Policy::Material;
            using Kernel = typename Policy::Kernel;
            using DofTraits = typename Policy::DofTraits;
            using DensityField = field::DensityField<Grid::Dim>;
            using Operator = fem::MatrixFreeOperator<Grid, DofTraits>;
            using DofMap = fem::DofMap<Grid, DofTraits>;
            using Preconditioner = JacobiPreconditioner<Operator>;

            using MaterialTensor = typename Material::MaterialTensor;

            /// @brief Element stiffness matrix type.
            using ElementStiffnessMatrix = typename Kernel::StiffnessMatrix;

            /// @brief Element field matrix type.
            using ElementFieldMatrix = typename Kernel::FieldMatrix;

            /// @brief Global field matrix type.
            using FieldMatrix = typename Policy::DofFieldMatrix;

            /// @brief Homogenized results type.
            using Results = typename Policy::Results;

            /**
             * @brief Solves the periodic cell problem.
             *
             * The homogenized material tensor M̄ is then computed via the physics-specific
             * homogenization functional:
             *
             * M̄=1/V∑ₑXₑᵀKₑXₑ
             *
             * @param[in] grid Grid.
             * @param[in] densityField Per-element density field defined on `grid`.
             * @param[in] material Base material.
             * @param[in] options Solver options.
             *
             * @returns Homogenized material tensor.
             *
             * @throws std::invalid_argument if `grid` and `densityField` do not
             * have the same resolution.
             * @throws std::runtime_error if the iterative linear solver fails.
             */
            Results solve(const Grid &grid, const DensityField &densityField, const Material &material, const SolverOptions &options = {}) const {
                if (grid.resolution() != densityField.resolution()) {
                    throw std::invalid_argument("Grid resolution must match density field resolution.");
                }

                const auto referenceNodes = grid.elementNodes(0);
                const ElementStiffnessMatrix elementKReference = Kernel::lhs(material, referenceNodes);
                const ElementFieldMatrix elementFReference = Kernel::rhs(material, referenceNodes);

                const Operator K(grid, densityField, elementKReference);
                const DofMap dofMap = K.dofMap();

                const FieldMatrix reducedF = buildReducedRhs_(grid, densityField, dofMap, elementFReference);
                const FieldMatrix reducedXMicro = solveReducedSystem_(K, reducedF, options);

                const FieldMatrix XMicro = expandReducedLhs_(grid, dofMap, reducedXMicro);
                const FieldMatrix XMacro = Policy::makeMacroscopicFields(grid);
                const FieldMatrix X = XMacro + XMicro;

                const MaterialTensor MBar = homogenize_(grid, densityField, elementKReference, X);

                return Policy::makeResults(MBar, X, XMacro, XMicro, options.fields);
            }

        private:
            /**
             * @brief Builds the reduced source matrix for the microscopic solve.
             *
             * Constructs the right-hand side matrix F in
             *
             * ```text
             * KX̃=F
             * ```
             *
             * on the reduced periodic dof space.
             *
             * @param[in] grid Grid.
             * @param[in] densityField Per-element density field defined on `grid`.
             * @param[in] dofMap Reduced periodic dof map.
             * @param[in] elementFReference Reference element source matrix for unit density.
             *
             * @returns Reduced periodic source matrix.
             */
            FieldMatrix buildReducedRhs_(const Grid &grid, const DensityField &densityField, const DofMap &dofMap, const ElementFieldMatrix &elementFReference) const noexcept {
                const int numReducedDofs = static_cast<int>(dofMap.numReducedDofs());
                const int numLoadCases = static_cast<int>(Policy::NumLoadCases);

                FieldMatrix reducedF = FieldMatrix::Zero(numReducedDofs, numLoadCases);

                for (std::size_t i = 0; i < grid.numElements(); ++i) {
                    const auto &reducedDofs = dofMap.reducedDofs(i);
                    const double density = densityField.getDensity(i);
                    const ElementFieldMatrix elementF = density * elementFReference;

                    for (std::size_t j = 0; j < reducedDofs.size(); ++j) {
                        const int reducedDof = reducedDofs[j];
                        const int localReducedDof = static_cast<int>(j);

                        if (reducedDof >= 0) {
                            reducedF.row(reducedDof) += elementF.row(localReducedDof);
                        }
                    }
                }

                return reducedF;
            }

            /**
             * @brief Solves the reduced microscopic system.
             *
             * Solves
             *
             * ```text
             * KX̃=F
             * ```
             *
             * on the reduced periodic dof space.
             *
             * @param[in] K Matrix-free operator for the reduced periodic stiffness matrix.
             * @param[in] reducedF Reduced periodic source matrix.
             * @param[in] options Solver options.
             *
             * @returns Reduced microscopic field matrix.
             *
             * @throws std::runtime_error if the iterative solver fails.
             */
            FieldMatrix solveReducedSystem_(const Operator &K, const FieldMatrix &reducedF, const SolverOptions &options) const {
                Eigen::ConjugateGradient<Operator, Eigen::Lower | Eigen::Upper, Preconditioner> cg;

                cg.setMaxIterations(options.maxIterations);
                cg.setTolerance(options.tolerance);
                cg.compute(K);

                FieldMatrix reducedXMicro = cg.solve(reducedF);

                if (cg.info() != Eigen::Success) {
                    std::string message = "Iterative solver failed";

                    if (cg.info() == Eigen::NoConvergence) {
                        message += ": No convergence before the iteration limit.";
                    }
                    else if (cg.info() == Eigen::NumericalIssue) {
                        message += ": Numerical issue encountered during the solve.";
                    }

                    throw std::runtime_error(message);
                }

                return reducedXMicro;
            }

            /**
             * @brief Expands the reduced periodic microscopic field matrix to the global dof space.
             *
             * Expands the left-hand side matrix X̃ in
             *
             * ```text
             * KX̃=F
             * ```
             *
             * to the global dof space.
             *
             * @param[in] grid Grid.
             * @param[in] dofMap Reduced periodic dof map.
             * @param[in] reducedX Reduced microscopic field matrix.
             *
             * @returns Global microscopic field matrix.
             */
            FieldMatrix expandReducedLhs_(const Grid &grid, const DofMap &dofMap, const FieldMatrix &reducedX) const noexcept {
                // Expand reduced periodic dofs to periodic dofs
                const std::size_t numPeriodicNodes = grid.numPeriodicNodes();
                const std::size_t numPeriodicDofs = DofTraits::NumNodeDofs * numPeriodicNodes;

                FieldMatrix periodicX = FieldMatrix::Zero(static_cast<int>(numPeriodicDofs), static_cast<int>(Policy::NumLoadCases));

                for (std::size_t i = 0; i < dofMap.numReducedDofs(); ++i) {
                    const std::size_t reducedDof = i;
                    const std::size_t periodicDof = dofMap.reducedToPeriodicDof(i);

                    periodicX.row(static_cast<int>(periodicDof)) = reducedX.row(static_cast<int>(reducedDof));
                }

                // Expand periodic dofs to global dofs
                const std::size_t numNodes = grid.numNodes();
                const std::size_t numDofs = DofTraits::NumNodeDofs * numNodes;

                FieldMatrix X = FieldMatrix::Zero(static_cast<int>(numDofs), static_cast<int>(Policy::NumLoadCases));

                for (std::size_t i = 0; i < grid.numElements(); ++i) {
                    const auto element = grid.element(i);
                    const auto periodicElement = grid.periodicElement(i);

                    const auto globalDofs = DofTraits::elementDofs(element, numNodes);
                    const auto periodicDofs = DofTraits::elementDofs(periodicElement, numPeriodicNodes);

                    for (std::size_t j = 0; j < globalDofs.size(); ++j) {
                        const std::size_t globalDof = globalDofs[j];
                        const std::size_t periodicDof = periodicDofs[j];

                        X.row(static_cast<int>(globalDof)) = periodicX.row(static_cast<int>(periodicDof));
                    }
                }

                return X;
            }

            /**
             * @brief Homogenized material tensor.
             *
             * Using the total field matrix X, the homogenized material tensor M̄
             * is computed  via the Hill-Mandel lemma:
             *
             * ```text
             * M̄=1/V∑ₑXₑᵀKₑXₑ
             * ```
             *
             * @param[in] grid Grid.
             * @param[in] densityField Per-element density field defined on `grid`.
             * @param[in] elementKReference Reference element stiffness matrix for unit density.
             * @param[in] X Global total field matrix.
             *
             * @returns Homogenized material tensor.
             */
            MaterialTensor homogenize_(const Grid &grid, const DensityField &densityField, const ElementStiffnessMatrix &elementKReference, const FieldMatrix &X) const {
                MaterialTensor MBar = MaterialTensor::Zero();

                const std::size_t numNodes = grid.numNodes();

                for (std::size_t i = 0; i < grid.numElements(); ++i) {
                    const double density = densityField.getDensity(i);
                    const auto elementK = density * elementKReference;
                    const auto element = grid.element(i);
                    const auto dofs = DofTraits::elementDofs(element, numNodes);

                    ElementFieldMatrix elementX;

                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const std::size_t dof = dofs[j];
                        const std::size_t localDof = j;
                        elementX.row(static_cast<int>(localDof)) = X.row(static_cast<int>(dof));
                    }

                    MBar.noalias() += elementX.transpose() * elementK * elementX;
                }

                MBar /= grid.measure();

                // Remove numerical asymmetry
                detail::symmetrize(MBar);

                return MBar;
            }
        };

    } // namespace solver

} // namespace monad
