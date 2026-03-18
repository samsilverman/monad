#pragma once

#include <cstddef>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "monad/field/density_field.hpp"
#include "monad/fem/operator/dof_map.hpp"
#include "monad/detail/eigen_utils.hpp"

namespace monad {

    namespace fem {

        /**
         * @brief Matrix-free operator for the reduced global stiffness matrix.
         *
         * This class applies the global stiffness operator without explicitly
         * assembling the global matrix.
         *
         * @tparam Grid Grid type (e.g. Quad4Grid).
         * @tparam Traits Dof traits type (e.g. LinearElasticDofTraits<2>).
         *
         * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
         */
        template <class Grid, class Traits>
        class MatrixFreeOperator : public Eigen::EigenBase<MatrixFreeOperator<Grid, Traits>> {
        public:
            using Element = typename Grid::Element;

            /// @brief Number of dofs per element.
            static constexpr int NumElementDofs = Element::NumNodes * Traits::NumNodeDofs;

            /// @brief Element dof vector type.
            using DofVector = Eigen::Vector<double, NumElementDofs>;

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
             * @brief Constructs a matrix-free operator.
             *
             * @param[in] grid Grid.
             * @param[in] densityField Per-element density field defined on `grid`.
             * @param[in] elementKReference Reference element stiffness matrix for unit density.
             *
             * @throws std::invalid_argument if `grid` and `densityField` do not have the same resolution.
             */
            MatrixFreeOperator(const Grid &grid, const field::DensityField<Grid::Dim> &densityField, const StiffnessMatrix &elementKReference)
                : densityField_(densityField), dofMap_(grid), elementKReference_(elementKReference) {
                if (grid.resolution() != densityField_.resolution()) {
                    throw std::invalid_argument( "Grid resolution must match density field resolution.");
                }
            }

            /// @brief Number of rows.
            int rows() const noexcept {
                return static_cast<int>(dofMap_.numReducedDofs());
            }

            /// @brief Number of columns.
            int cols() const noexcept {
                return static_cast<int>(dofMap_.numReducedDofs());
            }

            /// @brief `true` if the matrix-free operator is symmetric, `false` otherwise.
            bool isSymmetric() const noexcept {
                return detail::isSymmetric(elementKReference_);
            }

            /// @brief `true` if the matrix-free operator is positive semi-definite, `false` otherwise.
            bool isPSD() const noexcept {
                return detail::isPSD(elementKReference_);
            }

            /// @brief Per-element density field.
            const field::DensityField<Grid::Dim> &densityField() const noexcept {
                return densityField_;
            };

            /// @brief Reduced periodic dof map.
            const DofMap<Grid, Traits> &dofMap() const noexcept {
                return dofMap_;
            }

            /// @brief Reference element stiffness matrix for unit density.
            const StiffnessMatrix& elementKReference() const noexcept {
                return elementKReference_;
            }

            /**
             * @brief Matrix-vector multiplication Kx=y.
             *
             * @tparam Lhs Type of the left-hand side vector x.
             * @tparam Rhs Type of the right-hand side vector y.
             *
             * @param[in] lhs Left-hand side vector x.
             * @param[in,out] rhs Right-hand side vector y.
             */
            template <typename Lhs, typename Rhs>
            void apply(const Eigen::MatrixBase<Lhs> &lhs, Rhs &rhs) const {
#ifdef _OPENMP
#ifdef DEFAULT_OMP_NUM_THREADS
                omp_set_num_threads(DEFAULT_OMP_NUM_THREADS);
#endif
                #pragma omp parallel for schedule(static)
#endif
                for (std::size_t i = 0; i < dofMap_.size(); ++i) {
                    const auto &reducedDofs = dofMap_.reducedDofs(i);
                    const double density = densityField_.getDensity(i);

                    // Gather the element dof vector from the reduced periodic vector.
                    // Fixed dofs are represented by -1 in the map and are treated as zero.
                    DofVector x;

                    for (std::size_t j = 0; j < reducedDofs.size(); ++j) {
                        const int reducedDof = reducedDofs[j];
                        const int localReducedDof = static_cast<int>(j);

                        x(localReducedDof) = (reducedDof >= 0) ? lhs(reducedDof) : 0.0;
                    }

                    // Apply the density-scaled reference element stiffness
                    const DofVector y = density * (elementKReference_ * x);

                    // Scatter the local contribution back into the reduced global vector
                    for (std::size_t j = 0; j < reducedDofs.size(); ++j) {
                        const int reducedDof = reducedDofs[j];
                        const int localReducedDofs = static_cast<int>(j);

                        if (reducedDof >= 0) {
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            rhs(reducedDof) += y(localReducedDofs);
                        }
                    }
                }
            }

            /**
             * @brief Overloaded operator for matrix-vector multiplication y=Kx.
             *
             * @tparam Lhs Type of the right-hand side vector x.
             *
             * @param[in] x Left-hand side vector x.
             *
             * @returns Matrix-vector product Kx.
             */
            template <typename Lhs>
            Eigen::Product<MatrixFreeOperator, Lhs, Eigen::AliasFreeProduct>
            operator*(const Eigen::MatrixBase<Lhs> &x) const {
                return Eigen::Product<MatrixFreeOperator, Lhs, Eigen::AliasFreeProduct>(*this, x.derived());
            }

        private:
            /// @brief Per-element density field.
            const field::DensityField<Grid::Dim> &densityField_;

            /// @brief Reduced periodic dof map.
            DofMap<Grid, Traits> dofMap_;

            /// @brief Reference element stiffness matrix for unit density.
            const StiffnessMatrix &elementKReference_;
        };

    } // namespace fem

} // namespace monad

namespace Eigen {

    namespace internal {

        /**
         * @brief Eigen traits specialization for `MatrixFreeOperator`.
         *
         * This tells Eigen to treat `MatrixFreeOperator` like a sparse matrix in
         * expression templates.
         *
         * @tparam Grid Grid type (e.g. Quad4Grid).
         * @tparam Traits Dof traits type (e.g. LinearElasticDofTraits<2>).
         *
         * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
         */
        template <class Grid, class Traits>
        struct traits<monad::fem::MatrixFreeOperator<Grid, Traits>> : public Eigen::internal::traits<Eigen::SparseMatrix<double>> {};

        /**
         * @brief  Eigen product specialization for matrix-vector multiplication.
         *
         * @tparam Grid Grid type (e.g. Quad4Grid).
         * @tparam Traits Dof traits type (e.g. LinearElasticDofTraits<2>).
         * @tparam Lhs Type of the left-hand side vector.
         *
         * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
         */
        template <class Grid, class Traits, typename Lhs>
        struct generic_product_impl<
            monad::fem::MatrixFreeOperator<Grid, Traits>,
            Lhs,
            SparseShape,
            DenseShape,
            GemvProduct // General matrix-vector product
            > : generic_product_impl_base<monad::fem::MatrixFreeOperator<Grid, Traits>, Lhs,
                                          generic_product_impl<monad::fem::MatrixFreeOperator<Grid, Traits>, Lhs>> {
            using Operator = typename monad::fem::MatrixFreeOperator<Grid, Traits>;
            using Scalar = typename Product<Operator, Lhs>::Scalar;

            /**
             * @brief Matrix-vector multiplication Kx=y.
             *
             * @tparam Rhs Type of the right-hand side vector y.
             *
             * @param[in,out] rhs Right-hand side vector y.
             * @param[in] K Matrix-free operator.
             * @param[in] lhs Left-hand side vector x.
             * @param[in] alpha Scaling factor (must be 1).
             *
             * @note Inspiration: https://libeigen.gitlab.io/eigen/docs-nightly/group__MatrixfreeSolverExample.html
             */
            template <typename Rhs>
            static void scaleAndAddTo(Rhs &rhs, const Operator &K, const Lhs &lhs, const Scalar &alpha) {
                eigen_assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                K.apply(lhs, rhs);
            }
        };

    } // namespace internal

} // namespace Eigen
