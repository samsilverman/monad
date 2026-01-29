#include <numeric>
#include <stdexcept>
#include "monad/detail/mean.hpp"

namespace monad {

    namespace detail {

        double arithmeticMean(const std::vector<double> &x) noexcept {
            const double sum = std::accumulate(x.begin(), x.end(), 0.0);
            const double n = static_cast<double>(x.size());

            return sum / n;
        }

        double harmonicMean(const std::vector<double> &x) {
            const double invSum = std::accumulate(x.begin(), x.end(), 0.0, [](double accumulator, double value) {
                if (value == 0.0) {
                    throw std::invalid_argument("Scalar values must be nonzero.");
                }
                return accumulator + (1.0 / value);
            });
            const double n = static_cast<double>(x.size());

            return n / invSum;
        }

    } // namepsace detail

} // namespace monad
