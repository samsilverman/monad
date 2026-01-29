#pragma once

#include <vector>

namespace monad {

    namespace detail {

        /**
         * @brief Arithmetic mean for a list of scalar values.
         *
         * Given a collection of values {x₁,x₂,...,xₙ}, the arithmetic mean is defined as:
         *
         * (1/n)∑ᵢxᵢ
         *
         * @param[in] x List of scalar values.
         *
         * @returns Arithmetic mean for a list of scalar values.
         */
        double arithmeticMean(const std::vector<double> &x) noexcept;

        /**
         * @brief Harmonic mean for a list of scalar values.
         *
         * Given a collection of values {x₁,x₂,...,xₙ}, the harmonic mean is defined as:
         *
         * n/∑ᵢ(1/xᵢ)
         *
         * @param[in] x List of scalar values.
         *
         * @returns Harmonic mean for a list of scalar values.
         *
         * @throws std::invalid_argument if any entry in `x` is zero.
         */
        double harmonicMean(const std::vector<double> &x);

    } // namepsace detail

} // namespace monad
