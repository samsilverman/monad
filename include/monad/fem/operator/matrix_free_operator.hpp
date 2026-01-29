#pragma once

#include <cstddef>
#include <array>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "monad/detail/constants.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/grid/grid_base.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief Matrix-free operator for the global stiffness matrix.
         *
         * The operator represents the action of the global stiffness matrix
         * on a vector without explicitly assembling the matrix and is defined
         * on the reduced set of unconstrained periodic dofs.
         *
         * @tparam Grid Grid class (e.g. Quad4Grid).
         * @tparam Element Element class (e.g. Quad4).
         * @tparam Traits Matrix free operator traits class (e.g. LinearElasticMatrixFreeOperatorTraits<2>).
         *
         * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
         */
        template <class Grid, class Element, class Traits>
        class MatrixFreeOperator : public Eigen::EigenBase<MatrixFreeOperator<Grid, Element, Traits>> {
        public:
            static_assert(Element::Dim == Grid::Dim, "Element spatial dimension must equal grid spatial dimension.");

            /// @brief Number of dofs per element.
            static constexpr int NumElementDofs = Element::NumNodes * Traits::NumNodeDofs;

            using ElementsList = typename GridBase<Grid, Element>::ElementsList;
            using DensityList = typename GridBase<Grid, Element>::DensityList;
            using ElementDofList = std::vector<std::array<int, NumElementDofs>>;

            /// @brief Element stiffness matrix type.
            using StiffnessMatrix = Eigen::Matrix<double, NumElementDofs, NumElementDofs>;

            typedef double Scalar;
            typedef double RealScalar;
            typedef int StorageIndex;
            enum {
                ColsAtCompileTime = Eigen::Dynamic,
                MaxColsAtCompileTime = Eigen::Dynamic,
                IsRowMajor = false
            };

            /**
             * @brief Constructs the matrix-free operator.
             *
             * @param[in] grid Periodic unit cell grid.
             * @param[in] elementKReference Reference element stiffness matrix for unit density.
             */
            MatrixFreeOperator(const GridBase<Grid, Element> &grid, const StiffnessMatrix &elementKReference)
                : elementKReference_(elementKReference), elements_(grid.periodicElements()), densities_(grid.densities()) {
                const std::size_t numPeriodicNodes = grid.numPeriodicNodes();
                const std::size_t numReducedDofs = Traits::NumNodeDofs * numPeriodicNodes - Traits::NumFixedDofs;

                numRows_ = static_cast<int>(numReducedDofs);
                numCols_ = numRows_;

                // Precompute global-reduced dof mapping
                // Fixed dofs are marked with -1 and skipped during gather/scatter
                const std::size_t numElements = elements_.size();
                elementDofs_.resize(numElements);

                for (std::size_t i = 0; i < numElements; ++i) {
                    auto &element = elements_[i];
                    const auto dofs = Traits::dofs(element, numPeriodicNodes);

                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const std::size_t dof = dofs[j];
                        if (Traits::isFixedDof(dof, numPeriodicNodes)) {
                            elementDofs_[i][j] = -1;
                        }
                        else {
                            const std::size_t reducedDof = Traits::reducedDof(dof, numPeriodicNodes);
                            elementDofs_[i][j] = static_cast<int>(reducedDof);
                        }
                    }
                }
            }

            /// @brief Reference element stiffness matrix for unit density.
            const StiffnessMatrix &elementKReference() const noexcept {
                return elementKReference_;
            }

            /// @brief Periodic node indices for all elements.
            const ElementsList &elements() const noexcept {
                return elements_;
            }

            /// @brief Material densities for all elements.
            const DensityList &densities() const noexcept {
                return densities_;
            }

            /**
             * @brief Reduced dofs for all elements.
             *
             * @note Fixed dofs are set to -1.
             */
            const ElementDofList &elementDofs() const noexcept {
                return elementDofs_;
            }

            /// @brief Number of rows.
            int rows() const noexcept {
                return numRows_;
            }

            /// @brief Number of columns.
            int cols() const noexcept {
                return numCols_;
            }

            /// @brief `true` if the matrix-free operator is symmetric, `false` otherwise.
            bool isSymmetric() const noexcept {
                return detail::isSymmetric(elementKReference_);
            }

            /// @brief `true` if the matrix-free operator is positive semi-definite, `false` otherwise.
            bool isPSD() const noexcept {
                return detail::isPSD(elementKReference_);
            }

            /**
             * @brief Overloads the multiplication operator to define the matrix-vector product y=Kx.
             *
             * @tparam Rhs Type of the right-hand side vector.
             *
             * @param[in] x Right-hand side vector.
             *
             * @returns Product Kx.
             */
            template <typename Rhs>
            Eigen::Product<MatrixFreeOperator, Rhs, Eigen::AliasFreeProduct>
            operator*(const Eigen::MatrixBase<Rhs> &x) const {
                return Eigen::Product<MatrixFreeOperator, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
            }

        private:
            /// @brief Reference element stiffness matrix for unit density.
            const StiffnessMatrix &elementKReference_;

            /// @brief Number of rows.
            int numRows_;

            /// @brief Number of columns.
            int numCols_;

            /// @brief Periodic node indices for all elements.
            ElementsList elements_;

            /// @brief Material densities for all elements.
            DensityList densities_;

            /**
             * @brief Reduced dofs for all elements.
             *
             * @note Fixed dofs are set to -1.
             */
            ElementDofList elementDofs_;
        };

    } // namespace fem

} // namespace monad

namespace Eigen {

    namespace internal {

        /**
         * @brief This step informs Eigen's expression system that MatrixFreeOperator should be treated similarly to Eigen::SparseMatrix.
         *
         * @tparam Grid Grid class (e.g. Quad4Grid).
         * @tparam Element Element class (e.g. Quad4).
         * @tparam Traits Matrix free operator traits class (e.g. LinearElasticMatrixFreeOperatorTraits).
         *
         * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
         */
        template <class Grid, class Element, class Traits>
        struct traits<monad::fem::MatrixFreeOperator<Grid, Element, Traits>> : public Eigen::internal::traits<Eigen::SparseMatrix<double>> {};

        /**
         * @brief Specializes Eigen's product implementation for the matrix-vector multiplication Kx.
         *
         * @tparam Grid Grid class (e.g. Quad4Grid).
         * @tparam Element Element class (e.g. Quad4).
         * @tparam Traits Matrix free operator traits class (e.g. LinearElasticMatrixFreeOperatorTraits).
         * @tparam Rhs Type of the right-hand side vector.
         *
         * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
         */
        template <class Grid, class Element, class Traits, typename Rhs>
        struct generic_product_impl<
            monad::fem::MatrixFreeOperator<Grid, Element, Traits>,
            Rhs,
            SparseShape,
            DenseShape,
            GemvProduct // General matrix-vector product
            > : generic_product_impl_base<monad::fem::MatrixFreeOperator<Grid, Element, Traits>, Rhs,
                                          generic_product_impl<monad::fem::MatrixFreeOperator<Grid, Element, Traits>, Rhs>> {
            using Operator = typename monad::fem::MatrixFreeOperator<Grid, Element, Traits>;
            using Scalar = typename Product<Operator, Rhs>::Scalar;

            /**
             * @brief Implements the core matrix-vector multiplication y=Kx.
             *
             * @tparam Dest The type of the destination vector y.
             *
             * @param[in,out] dst Destination vector y.
             * @param[in] lhs Matrix-free operator K.
             * @param[in] rhs Right-hand side vector x.
             * @param[in] alpha Scaling factor (must be 1).
             *
             * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
             */
            template <typename Dest>
            static void scaleAndAddTo(Dest &dst, const Operator &lhs, const Rhs &rhs, const Scalar &alpha) {
                eigen_assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                using DofVector = Eigen::Vector<double, Operator::NumElementDofs>;

                const auto &elementKReference = lhs.elementKReference();
                const auto &elements = lhs.elements();
                const auto &densities = lhs.densities();
                const auto &elementDofs = lhs.elementDofs();

#ifdef _OPENMP
#ifdef DEFAULT_OMP_NUM_THREADS
                omp_set_num_threads(DEFAULT_OMP_NUM_THREADS);
#endif
                #pragma omp parallel for schedule(static)
#endif
                for (std::size_t i = 0; i < elements.size(); ++i) {
                    const double density = densities[i];
                    const auto &dofs = elementDofs[i];

                    // --- Gather ---
                    // Assemble the local vector x from the reduced global vector rhs.
                    // Fixed dofs are treated as zero.
                    DofVector x;

                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const int globalDof = dofs[j];
                        const int localDof = static_cast<int>(j);

                        if (globalDof >= 0) {
                            x(localDof) = rhs(globalDof);
                        }
                        else {
                            x(localDof) = 0.0;
                        }
                    }

                    // --- Apply ---
                    // Apply the element stiffness Kx (scaled by material density).
                    const DofVector y = density * (elementKReference * x);

                    // --- Scatter ---
                    // Accumulate local element Kx contributions into the reduced global vector dst.
                    for (std::size_t j = 0; j < dofs.size(); ++j) {
                        const int globalDof = dofs[j];
                        const int localDof = static_cast<int>(j);

                        if (globalDof >= 0) {

#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            dst(globalDof) += y(localDof);
                        }
                    }
                }
            }
        };

    } // namespace internal

} // namespace Eigen
