#pragma once

#include <Eigen/Core>

namespace monad {

    namespace fem {

        namespace mechanical {

            /**
             * @brief Linear elastic kernel traits.
             *
             * This struct isolates dimension-specific logic so that the finite element kernel
             * can be written generically in terms of spatial dimension.
             *
             * @tparam D Spatial dimension (2 or 3).
             */
            template <int D>
            struct LinearElasticKernelTraits;

            /**
             * @brief 2D linear elastic kernel traits.
             *
             * This struct isolates 2D-specific logic so that the finite element kernel
             * can be written generically in terms of spatial dimension.
             */
            template <>
            struct LinearElasticKernelTraits<2> {
                /**
                 * @brief 2D element B matrix evaluated at a local point.
                 *
                 * The B matrix relates an element's nodal displacements u to strains ε in Voigt
                 * notation at a local point:
                 *
                 * ε=Bu
                 *
                 * @tparam BMatrix B matrix type.
                 * @tparam ShapeFuncGradMatrix Shape function gradient matrix type.
                 *
                 * @param[in,out] B Element B matrix evaluated at a point.
                 * @param[in] dN Element shape function gradients evaluated at a point.
                 */
                template <class BMatrix, class ShapeFuncGradMatrix>
                static void fillB(BMatrix &B, const ShapeFuncGradMatrix &dN) {
                    B.setZero();

                    // Normal strain ε₁₁
                    B(0, Eigen::seq(0, Eigen::indexing::last, 2)) = dN.row(0);

                    // Normal strain ε₂₂
                    B(1, Eigen::seq(1, Eigen::indexing::last, 2)) = dN.row(1);

                    // Shear strain ε₁₂
                    B(2, Eigen::seq(0, Eigen::indexing::last, 2)) = dN.row(1);
                    B(2, Eigen::seq(1, Eigen::indexing::last, 2)) = dN.row(0);
                }
            };

            /**
             * @brief 3D linear elastic kernel traits.
             *
             * This struct isolates 3D-specific logic so that the finite element kernel
             * can be written generically in terms of spatial dimension.
             */
            template <>
            struct LinearElasticKernelTraits<3> {
                /**
                 * @brief 3D element B matrix evaluated at a local point.
                 *
                 * The B matrix relates an element's nodal displacements u to strains ε in Voigt
                 * notation at a local point:
                 *
                 * ε=Bu
                 *
                 * @tparam BMatrix B matrix type.
                 * @tparam ShapeFuncGradMatrix Shape function gradient matrix type.
                 *
                 * @param[in,out] B Element B matrix evaluated at a point.
                 * @param[in] dN Element shape function gradients evaluated at a point.
                 */
                template <class BMatrix, class ShapeFuncGradMatrix>
                static void fillB(BMatrix &B, const ShapeFuncGradMatrix &dN) {
                    B.setZero();

                    // Normal strain ε₁₁
                    B(0, Eigen::seq(0, Eigen::indexing::last, 3)) = dN.row(0);

                    // Normal strain ε₂₂
                    B(1, Eigen::seq(1, Eigen::indexing::last, 3)) = dN.row(1);

                    // Normal strain ε₃₃
                    B(2, Eigen::seq(2, Eigen::indexing::last, 3)) = dN.row(2);

                    // Shear strain ε₁₂
                    B(3, Eigen::seq(0, Eigen::indexing::last, 3)) = dN.row(1);
                    B(3, Eigen::seq(1, Eigen::indexing::last, 3)) = dN.row(0);

                    // Shear strain ε₁₃
                    B(4, Eigen::seq(0, Eigen::indexing::last, 3)) = dN.row(2);
                    B(4, Eigen::seq(2, Eigen::indexing::last, 3)) = dN.row(0);

                    // Shear strain ε₂₃
                    B(5, Eigen::seq(1, Eigen::indexing::last, 3)) = dN.row(2);
                    B(5, Eigen::seq(2, Eigen::indexing::last, 3)) = dN.row(1);
                }
            };

        } // namespace mechanical

    } // namespace fem

} // namespace monad
