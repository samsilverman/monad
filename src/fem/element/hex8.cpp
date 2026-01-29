#include <cmath>
#include <algorithm>
#include "monad/fem/element/hex8.hpp"

namespace monad {

	namespace fem {

		Hex8::NodesMatrix Hex8::localNodes() noexcept {
			NodesMatrix nodes {
				{-1.0, -1.0, -1.0},
				{1.0, -1.0, -1.0},
				{1.0, 1.0, -1.0},
				{-1.0, 1.0, -1.0},
				{-1.0, -1.0, 1.0},
				{1.0, -1.0, 1.0},
				{1.0, 1.0, 1.0},
				{-1.0, 1.0, 1.0}
			};

			return nodes;
		}

		Hex8::ShapeFuncVector Hex8::shapeFunctions(const Point &point) noexcept {
			const double xi = point(0);
			const double eta = point(1);
			const double zeta = point(2);

			ShapeFuncVector N;
			N(0) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
			N(1) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta);
			N(2) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta);
			N(3) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta);
			N(4) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta);
			N(5) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta);
			N(6) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta);
			N(7) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta);

			return N;
		}

		Hex8::ShapeFuncGradMatrix Hex8::gradShapeFunctions(const Point &point) noexcept {
			const double xi = point(0);
			const double eta = point(1);
			const double zeta = point(2);

			ShapeFuncGradMatrix dN;

			dN.col(0) << -0.125 * (eta - 1.0) * (zeta - 1),
						 -0.125 * (xi - 1.0) * (zeta - 1),
						 -0.125 * (xi - 1.0) * (eta - 1);
			dN.col(1) << 0.125 * (eta - 1.0) * (zeta - 1),
						 0.125 * (xi + 1.0) * (zeta - 1),
						 0.125 * (xi + 1.0) * (eta - 1);
			dN.col(2) << -0.125 * (eta + 1.0) * (zeta - 1),
						 -0.125 * (xi + 1.0) * (zeta - 1),
						 -0.125 * (xi + 1.0) * (eta + 1);
			dN.col(3) << 0.125 * (eta + 1.0) * (zeta - 1),
						 0.125 * (xi - 1.0) * (zeta - 1),
						 0.125 * (xi - 1.0) * (eta + 1);
			dN.col(4) << 0.125 * (eta - 1.0) * (zeta + 1),
						 0.125 * (xi - 1.0) * (zeta + 1),
						 0.125 * (xi - 1.0) * (eta - 1);
			dN.col(5) << -0.125 * (eta - 1.0) * (zeta + 1),
						 -0.125 * (xi + 1.0) * (zeta + 1),
						 -0.125 * (xi + 1.0) * (eta - 1);
			dN.col(6) << 0.125 * (eta + 1.0) * (zeta + 1),
						 0.125 * (xi + 1.0) * (zeta + 1),
						 0.125 * (xi + 1.0) * (eta + 1);
			dN.col(7) << -0.125 * (eta + 1.0) * (zeta + 1),
						 -0.125 * (xi - 1.0) * (zeta + 1),
						 -0.125 * (xi - 1.0) * (eta + 1);

			return dN;
		}

		Hex8::QuadratureRule Hex8::quadratureRule() noexcept {
			const double a = 1.0 / std::sqrt(3.0);

			Hex8::QuadratureRule rule;

			rule.points = {
				Point(-a, -a, -a),
				Point(-a, a, -a),
				Point(a, -a, -a),
				Point(a, a, -a),
				Point(-a, -a, a),
				Point(-a, a, a),
				Point(a, -a, a),
				Point(a, a, a)
			};

			std::fill(rule.weights.begin(), rule.weights.end(), 1.0);

			return rule;
		}

		int Hex8::gmshElementType() noexcept {
			return 5;
		}

		Hex8::NodeIndicesList Hex8::gmshNodeOrdering() noexcept {
			return {0, 1, 5, 4, 3, 2, 6, 7};
		}

	} // namespace fem

} // namespace monad
