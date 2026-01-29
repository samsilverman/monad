#pragma once

#include <cstddef>
#include <Eigen/Core>
#include "monad/grid/grid_base.hpp"
#include "monad/fem/operator/mechanical/linear_elastic_matrix_free_operator_traits.hpp"

namespace monad {

    namespace solver {

        namespace mechanical {

            /**
             * @brief Linear elastic physics traits.
             *
             * This struct isolates dimension-specific logic so that linear elastic physics
             * can be written generically in terms of spatial dimension.
             *
             * @tparam D Spatial dimension (2 or 3).
             */
            template <int D>
            struct LinearElasticPhysicsTraits;

            /**
             * @brief 2D linear elastic physics traits.
             *
             * This struct isolates 2D-specific logic so that linear elastic physics
             * can be written generically in terms of spatial dimension.
             */
            template <>
            struct LinearElasticPhysicsTraits<2> {
                using OperatorTraits = fem::mechanical::LinearElasticMatrixFreeOperatorTraits<2>;

                /**
                 * @brief Macroscopic nodal displacements.
                 *
                 * The macroscopic displacements Ū are the displacements
                 * resulting from the macroscopic strains ε̄:
                 *
                 * Ū=ε̄x
                 *
                 * @tparam Grid Grid class (e.g. Quad4Grid).
                 * @tparam Element Element class (e.g. Quad4).
                 *
                 * @param[in] grid Periodic unit cell grid.
                 *
                 * @returns Macroscopic nodal displacements.
                 */
                template <class Grid, class Element>
                static Eigen::MatrixX3d macroDisplacements(const GridBase<Grid, Element> &grid) {
                    static_assert(Element::Dim == Grid::Dim, "Element spatial dimension must equal grid spatial dimension.");

                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = OperatorTraits::NumNodeDofs * numNodes;

                    Eigen::MatrixX3d U = Eigen::MatrixX3d::Zero(static_cast<int>(numDofs), 3);

                    for (std::size_t i = 0; i < grid.numNodes(); ++i) {
                        auto node = grid.node(i);

                        const double x = node(0);
                        const double y = node(1);

                        const int u = static_cast<int>(2 * i);
                        const int v = static_cast<int>(2 * i + 1);

                        U(u, 0) = x;
                        U(v, 1) = y;
                        U(u, 2) = 0.5 * y;
                        U(v, 2) = 0.5 * x;
                    }

                    return U;
                }
            };

            /**
             * @brief 3D linear elastic physics traits.
             *
             * This struct isolates 3D-specific logic so that linear elastic physics
             * can be written generically in terms of spatial dimension.
             */
            template <>
            struct LinearElasticPhysicsTraits<3> {
                using OperatorTraits = fem::mechanical::LinearElasticMatrixFreeOperatorTraits<3>;

                using MatrixX6d = Eigen::Matrix<double, Eigen::Dynamic, 6>;

                /**
                 * @brief Macroscopic nodal displacements.
                 *
                 * The macroscopic displacements Ū are the displacements
                 * resulting from the macroscopic strains ε̄:
                 *
                 * Ū=ε̄x
                 *
                 * @tparam Grid Grid class (e.g. Hex8Grid).
                 * @tparam Element Element class (e.g. Hex8).
                 *
                 * @param[in] grid Periodic unit cell grid.
                 *
                 * @returns Macroscopic nodal displacements.
                 */
                template <class Grid, class Element>
                static MatrixX6d macroDisplacements(const GridBase<Grid, Element> &grid) {
                    static_assert(Element::Dim == Grid::Dim, "Element spatial dimension must equal grid spatial dimension.");

                    const std::size_t numNodes = grid.numNodes();
                    const std::size_t numDofs = OperatorTraits::NumNodeDofs * numNodes;

                    MatrixX6d U = MatrixX6d::Zero(static_cast<int>(numDofs), 6);

                    for (std::size_t i = 0; i < grid.numNodes(); ++i) {
                        const auto node = grid.node(i);

                        const double x = node(0);
                        const double y = node(1);
                        const double z = node(2);

                        const int u = static_cast<int>(3 * i);
                        const int v = static_cast<int>(3 * i + 1);
                        const int w = static_cast<int>(3 * i + 2);

                        U(u, 0) = x;
                        U(v, 1) = y;
                        U(w, 2) = z;
                        U(u, 3) = 0.5 * y;
                        U(v, 3) = 0.5 * x;
                        U(u, 4) = 0.5 * z;
                        U(w, 4) = 0.5 * x;
                        U(w, 5) = 0.5 * y;
                        U(v, 5) = 0.5 * z;
                    }

                    return U;
                }
            };

        } // namespace mechanical

    } // namespace solver

} // namespace monad
