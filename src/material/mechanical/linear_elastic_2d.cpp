#include <stdexcept>
#include <string>
#include <cmath>
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"

namespace monad {

    LinearElasticMaterial2d::LinearElasticMaterial2d(double E, double nu, PlaneCondition condition) {
        if (E <= 0.0) {
            throw std::invalid_argument("E (" + std::to_string(E) + ") must be positive.");
        }
        if (nu <= -1.0 || nu >= 0.5) {
            throw std::invalid_argument("nu (" + std::to_string(nu) + ") must be in range (-1,0.5).");
        }

        switch (condition) {
        case PlaneCondition::PlaneStress:
            C_ << 1.0, nu, 0.0,
                  nu, 1.0, 0.0,
                  0.0, 0.0, (1 - nu) / 2.0;
            C_ *= E / (1.0 - std::pow(nu, 2.0));
            break;
        case PlaneCondition::PlaneStrain:
            C_ << 1.0 - nu, nu, 0.0,
                  nu, 1.0 - nu, 0.0,
                  0.0, 0.0, (1.0 - 2.0 * nu) / 2.0;
            C_ *= E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            break;
        default:
            throw std::invalid_argument("Unknown plane condition.");
        }
    }

} // namespace monad
