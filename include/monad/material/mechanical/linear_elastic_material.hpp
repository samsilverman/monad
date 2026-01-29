#pragma once

#include <stdexcept>
#include <Eigen/Core>
#include "monad/detail/mean.hpp"
#include "monad/detail/eigen_utils.hpp"
#include "monad/grid/grid_base.hpp"

namespace monad {

    /**
     * @brief Represents a linear elastic material model.
     *
     * By Hooke's law, stress σ∈ℝᵛ is a linear function of strain ε∈ℝᵛ:
     *
     * σ=Cε
     *
     * - C∈Symᵥ(ℝ) is the stiffness tensor (in Voigt notation).
     *
     * @tparam D Spatial dimension (2 or 3).
     */
    template <int D>
    class LinearElasticMaterial {
    public:
        static_assert(D == 2 || D == 3, "Spatial dimension D must be 2 or 3.");

        /// @brief Spatial dimension (2 or 3).
        static constexpr int Dim = D;

        /// @brief Number of components in Voigt notation.
        static constexpr int VoigtSize = (Dim == 2) ? 3 : 6;

        /// @brief Stiffness tensor (in Voigt notation) type.
        using MaterialTensor = Eigen::Matrix<double, VoigtSize, VoigtSize>;

        /// @brief Default constructor.
        LinearElasticMaterial() = default;

        /**
         * @brief Constructs a linear elastic material.
         *
         * @param[in] C Stiffness tensor (in Voigt notation).
         *
         * @throws std::invalid_argument if `C` is not positive definite.
         */
        explicit LinearElasticMaterial(const MaterialTensor &C)
            : C_(C) {
            if (!detail::isPD(C_)) {
                throw std::invalid_argument("Stiffness tensor is not positive definite.");
            }
        }

        /// @brief Stiffness tensor C (in Voigt notation).
        const MaterialTensor &materialTensor() const noexcept {
            return C_;
        }

        /**
         * @brief Voigt upper bound for the homogenized stiffness tensor.
         *
         * @tparam Grid Grid class (e.g. Quad4Grid).
         * @tparam Element Element class (e.g. Quad4).
         *
         * @param[in] grid Periodic unit cell grid.
         *
         * @returns Voigt upper bound for the homogenized stiffness tensor.
         */
        template <class Grid, class Element>
        MaterialTensor voigt(const GridBase<Grid, Element> &grid) const noexcept {
            const double density = detail::arithmeticMean(grid.densities());

            return density * C_;
        }

        /**
         * @brief Reuss lower bound for the homogenized stiffness tensor.
         *
         * @tparam Grid Grid class (e.g. Quad4Grid).
         * @tparam Element Element class (e.g. Quad4).
         *
         * @param[in] grid Periodic unit cell grid.
         *
         * @returns Reuss lower bound for the homogenized stiffness tensor.
         */
        template <class Grid, class Element>
        MaterialTensor reuss(const GridBase<Grid, Element> &grid) const noexcept {
            const double density = detail::harmonicMean(grid.densities());

            return density * C_;
        }

        /// @brief Equality comparison.
        bool operator==(const LinearElasticMaterial &other) const noexcept {
            return C_.isApprox(other.C_);
        }

        /// @brief Inequality comparison.
        bool operator!=(const LinearElasticMaterial &other) const noexcept {
            return !(*this == other);
        }

    protected:
        /// @brief Stiffness tensor (in Voigt notation).
        MaterialTensor C_;
    };

} // namespace monad
