#include <cmath>
#include "monad/fem/element/quad8.hpp"

namespace monad {

	namespace fem {

		Quad8::NodesMatrix Quad8::localNodes() noexcept {
			NodesMatrix nodes {
				{-1.0, -1.0},
				{1.0, -1.0},
				{1.0, 1.0},
				{-1.0, 1.0},
				{0.0, -1.0},
				{1.0, 0.0},
				{0.0, 1.0},
				{-1.0, 0.0}
			};

			return nodes;
		}

		Quad8::ShapeFuncVector Quad8::shapeFunctions(const Point &point) noexcept {
			const double xi = point(0);
			const double eta = point(1);

			const double xi2 = std::pow(xi, 2);
			const double eta2 = std::pow(eta, 2);

			ShapeFuncVector N;
			N(0) = 0.25 * (1.0 - xi) * (1.0 - eta) * (-1.0 - xi - eta);
			N(1) = 0.25 * (1.0 + xi) * (1.0 - eta) * (-1.0 + xi - eta);
			N(2) = 0.25 * (1.0 + xi) * (1.0 + eta) * (-1.0 + xi + eta);
			N(3) = 0.25 * (1.0 - xi) * (1.0 + eta) * (-1.0 - xi + eta);
			N(4) = 0.5 * (1.0 - xi2) * (1.0 - eta);
			N(5) = 0.5 * (1.0 + xi) * (1.0 - eta2);
			N(6) = 0.5 * (1.0 - xi2) * (1.0 + eta);
			N(7) = 0.5 * (1.0 - xi) * (1.0 - eta2);

			return N;
		}

		Quad8::ShapeFuncGradMatrix Quad8::gradShapeFunctions(const Point &point) noexcept {
			const double xi = point(0);
			const double eta = point(1);

			const double xi2 = std::pow(xi, 2);
			const double eta2 = std::pow(eta, 2);

			ShapeFuncGradMatrix dN;

			dN.col(0) << -0.25 * (eta - 1.0) * (2.0 * xi + eta),
						 -0.25 * (xi - 1.0) * (xi + 2.0 * eta);
			dN.col(1) << -0.25 * (eta - 1.0) * (2.0 * xi - eta),
						 -0.25 * (xi + 1.0) * (xi - 2.0 * eta);
			dN.col(2) << 0.25 * (eta + 1.0) * (2.0 * xi + eta),
						 0.25 * (xi + 1.0) * (xi + 2.0 * eta);
			dN.col(3) << 0.25 * (eta + 1.0) * (2.0 * xi - eta),
						 0.25 * (xi - 1.0) * (xi - 2.0 * eta);
			dN.col(4) << xi * (eta - 1.0),
						 0.5 * (xi2 - 1.0);
			dN.col(5) << 0.5 * (1.0 - eta2),
						 -1.0 * (xi + 1.0) * eta;
			dN.col(6) << -1.0 * xi * (eta + 1.0),
						 0.5 * (1.0 - xi2);
			dN.col(7) << 0.5 * (eta2 - 1.0),
						 (xi - 1.0) * eta;

			return dN;
		}

		Quad8::QuadratureRule Quad8::quadratureRule() noexcept {
			const double a = std::sqrt(3.0 / 5.0);

			Quad8::QuadratureRule rule;

			rule.points = {
				Point(-a, -a),
				Point(-a, 0.0),
				Point(-a, a),
				Point(0.0, -a),
				Point(0.0, 0.0),
				Point(0.0, a),
				Point(a, -a),
				Point(a, 0.0),
				Point(a, a)
			};

			rule.weights = {
				25.0 / 81.0, // (5/9) * (5/9)
				40.0 / 81.0, // (5/9) * (8/9)
				25.0 / 81.0, // (5/9) * (5/9)
				40.0 / 81.0, // (8/9) * (5/9)
				64.0 / 81.0, // (8/9) * (8/9)
				40.0 / 81.0, // (8/9) * (5/9)
				25.0 / 81.0, // (5/9) * (5/9)
				40.0 / 81.0, // (5/9) * (8/9)
				25.0 / 81.0	 // (5/9) * (5/9)
			};

			return rule;
		}

		int Quad8::gmshElementType() noexcept {
			return 16;
		}

		Quad8::NodeIndicesList Quad8::gmshNodeOrdering() noexcept {
			return {0, 1, 2, 3, 4, 5, 6, 7};
		}

	} // namespace fem

} // namespace monad
