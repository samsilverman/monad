#pragma once

#include <vector>

namespace monad {

    namespace detail {

        /**
         * @brief Arithmetic mean of a list of scalar values.
         *
         * Given values
         *
         * ```text
         * {x₁, x₂, ..., xₙ}
         * ```
         * 
         * the arithmetic mean is defined as:
         *
         * ```text
         * (1/n)∑ᵢxᵢ
         * ```
         *
         * @param[in] x Values.
         *
         * @returns Arithmetic mean of `x`.
         */
        double arithmeticMean(const std::vector<double> &x) noexcept;

        /**
         * @brief Harmonic mean of a list of scalar values.
         *
         * Given values
         *
         * ```text
         * {x₁,x₂,...,xₙ}
         * ```
         *
         * the harmonic mean is defined as:
         *
         * ```text
         * n/∑ᵢ(1/xᵢ)
         * ```
         *
         * @param[in] x Values.
         *
         * @returns Harmonic mean of `x`.
         *
         * @throws std::invalid_argument if any entry of `x` is zero.
         */
        double harmonicMean(const std::vector<double> &x);

    } // namepsace detail

} // namespace monad
