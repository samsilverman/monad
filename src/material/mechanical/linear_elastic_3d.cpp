#include <stdexcept>
#include <string>
#include <cmath>
#include "monad/material/mechanical/linear_elastic_material_3d.hpp"

namespace monad {

    LinearElasticMaterial3d::LinearElasticMaterial3d(double E, double nu) {
        if (E <= 0.0) {
            throw std::invalid_argument("E (" + std::to_string(E) + ") must be positive.");
        }
        if (nu <= -1.0 || nu >= 0.5) {
            throw std::invalid_argument("nu (" + std::to_string(nu) + ") must be in range (-1,0.5).");
        }

        const double lame1 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double lame2 = E / (2.0 * (1.0 + nu));

        C_.setZero();

        C_(0, 0) = C_(1, 1) = C_(2, 2) = lame1 + 2.0 * lame2;

        C_(0, 1) = C_(1, 0) = lame1;
        C_(0, 2) = C_(2, 0) = lame1;
        C_(1, 2) = C_(2, 1) = lame1;

        C_(3, 3) = lame2;
        C_(4, 4) = lame2;
        C_(5, 5) = lame2;
    }

} // namespace monad
