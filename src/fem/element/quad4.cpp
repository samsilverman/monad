#include <cmath>
#include <algorithm>
#include "monad/fem/element/quad4.hpp"

namespace monad {

    namespace fem {

        Quad4::NodesMatrix Quad4::localNodes() noexcept {
            NodesMatrix nodes {
                {-1.0, -1.0},
                {1.0, -1.0},
                {1.0, 1.0},
                {-1.0, 1.0}
            };

            return nodes;
        }

        Quad4::ShapeFuncVector Quad4::shapeFunctions(const Point &point) noexcept {
            const double xi = point(0);
            const double eta = point(1);

            ShapeFuncVector N;
            N(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
            N(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
            N(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
            N(3) = 0.25 * (1.0 - xi) * (1.0 + eta);

            return N;
        }

        Quad4::ShapeFuncGradMatrix Quad4::gradShapeFunctions(const Point &point) noexcept {
            const double xi = point(0);
            const double eta = point(1);

            ShapeFuncGradMatrix dN;

            dN.col(0) << 0.25 * -(1.0 - eta),
                         0.25 * -(1.0 - xi);
            dN.col(1) << 0.25 * (1.0 - eta),
                         0.25 * -(1.0 + xi);
            dN.col(2) << 0.25 * (1.0 + eta),
                         0.25 * (1.0 + xi);
            dN.col(3) << 0.25 * -(1.0 + eta),
                         0.25 * (1.0 - xi);

            return dN;
        }

        Quad4::QuadratureRule Quad4::quadratureRule() noexcept {
            const double a = 1.0 / std::sqrt(3.0);

            Quad4::QuadratureRule rule;

            rule.points = {
                Point(-a, -a),
                Point(-a, a),
                Point(a, -a),
                Point(a, a)
            };

            std::fill(rule.weights.begin(), rule.weights.end(), 1.0);

            return rule;
        }

        int Quad4::gmshElementType() noexcept {
            return 3;
        }

        Quad4::NodeIndicesList Quad4::gmshNodeOrdering() noexcept {
            return {0, 1, 2, 3};
        }

    } // namespace fem

} // namespace monad
