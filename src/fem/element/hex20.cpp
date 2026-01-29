#include <cmath>
#include "monad/fem/element/hex20.hpp"

namespace monad {

	namespace fem {

		Hex20::NodesMatrix Hex20::localNodes() noexcept {
			NodesMatrix nodes {
				{-1.0, -1.0, -1.0},
				{1.0, -1.0, -1.0},
				{1.0, 1.0, -1.0},
				{-1.0, 1.0, -1.0},
				{-1.0, -1.0, 1.0},
				{1.0, -1.0, 1.0},
				{1.0, 1.0, 1.0},
				{-1.0, 1.0, 1.0},
				{0.0, -1.0, -1.0},
				{1.0, 0.0, -1.0},
				{0.0, 1.0, -1.0},
				{-1.0, 0.0, -1.0},
				{0.0, -1.0, 1.0},
				{1.0, 0.0, 1.0},
				{0.0, 1.0, 1.0},
				{-1.0, 0.0, 1.0},
				{-1.0, -1.0, 0.0},
				{1.0, -1.0, 0.0},
				{1.0, 1.0, 0.0},
				{-1.0, 1.0, 0.0}
			};

			return nodes;
		}

		Hex20::ShapeFuncVector Hex20::shapeFunctions(const Point &point) noexcept {
			const double xi = point(0);
			const double eta = point(1);
			const double zeta = point(2);

			const double xi2 = std::pow(xi, 2);
			const double eta2 = std::pow(eta, 2);
			const double zeta2 = std::pow(zeta, 2);

			ShapeFuncVector N;
			N(0) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1 - zeta) * (-2.0 - xi - eta - zeta);
			N(1) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1 - zeta) * (-2.0 + xi - eta - zeta);
			N(2) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1 - zeta) * (-2.0 + xi + eta - zeta);
			N(3) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1 - zeta) * (-2.0 - xi + eta - zeta);
			N(4) = 0.125 * (1.0 - xi) * (1.0 - eta) * (1 + zeta) * (-2.0 - xi - eta + zeta);
			N(5) = 0.125 * (1.0 + xi) * (1.0 - eta) * (1 + zeta) * (-2.0 + xi - eta + zeta);
			N(6) = 0.125 * (1.0 + xi) * (1.0 + eta) * (1 + zeta) * (-2.0 + xi + eta + zeta);
			N(7) = 0.125 * (1.0 - xi) * (1.0 + eta) * (1 + zeta) * (-2.0 - xi + eta + zeta);
			N(8) = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 - zeta);
			N(9) = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 - zeta);
			N(10) = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 - zeta);
			N(11) = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 - zeta);
			N(12) = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 + zeta);
			N(13) = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 + zeta);
			N(14) = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 + zeta);
			N(15) = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 + zeta);
			N(16) = 0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta2);
			N(17) = 0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta2);
			N(18) = 0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta2);
			N(19) = 0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta2);

			return N;
		}

		Hex20::ShapeFuncGradMatrix Hex20::gradShapeFunctions(const Point &point) noexcept {
			const double xi = point(0);
			const double eta = point(1);
			const double zeta = point(2);

			const double xi2 = std::pow(xi, 2);
			const double eta2 = std::pow(eta, 2);
			const double zeta2 = std::pow(zeta, 2);

			ShapeFuncGradMatrix dN;

			dN.col(0) << -0.125 * (eta - 1.0) * (zeta - 1.0) * (-2.0 * xi - eta - zeta - 1.0),
						 -0.125 * (xi - 1.0) * (zeta - 1.0) * (-xi - 2.0 * eta - zeta - 1.0),
						 -0.125 * (xi - 1.0) * (eta - 1.0) * (-xi - eta - 2.0 * zeta - 1.0);
			dN.col(1) << 0.125 * (eta - 1.0) * (zeta - 1.0) * (2.0 * xi - eta - zeta - 1.0),
						 0.125 * (xi + 1.0) * (zeta - 1.0) * (xi - 2.0 * eta - zeta - 1.0),
						 0.125 * (xi + 1.0) * (eta - 1.0) * (xi - eta - 2.0 * zeta - 1.0);
			dN.col(2) << -0.125 * (eta + 1.0) * (zeta - 1.0) * (2.0 * xi + eta - zeta - 1.0),
						 -0.125 * (xi + 1.0) * (zeta - 1.0) * (xi + 2.0 * eta - zeta - 1.0),
						 -0.125 * (xi + 1.0) * (eta + 1.0) * (xi + eta - 2.0 * zeta - 1.0);
			dN.col(3) << 0.125 * (eta + 1.0) * (zeta - 1.0) * (-2.0 * xi + eta - zeta - 1.0),
						 0.125 * (xi - 1.0) * (zeta - 1.0) * (-xi + 2.0 * eta - zeta - 1.0),
						 0.125 * (xi - 1.0) * (eta + 1.0) * (-xi + eta - 2.0 * zeta - 1.0);
			dN.col(4) << 0.125 * (eta - 1.0) * (zeta + 1.0) * (-2.0 * xi - eta + zeta - 1.0),
						 0.125 * (xi - 1.0) * (zeta + 1.0) * (-xi - 2.0 * eta + zeta - 1.0),
						 0.125 * (xi - 1.0) * (eta - 1.0) * (-xi - eta + 2.0 * zeta - 1.0);
			dN.col(5) << -0.125 * (eta - 1.0) * (zeta + 1.0) * (2.0 * xi - eta + zeta - 1.0),
						 -0.125 * (xi + 1.0) * (zeta + 1.0) * (xi - 2.0 * eta + zeta - 1.0),
						 -0.125 * (xi + 1.0) * (eta - 1.0) * (xi - eta + 2.0 * zeta - 1.0);
			dN.col(6) << 0.125 * (eta + 1.0) * (zeta + 1.0) * (2.0 * xi + eta + zeta - 1.0),
						 0.125 * (xi + 1.0) * (zeta + 1.0) * (xi + 2.0 * eta + zeta - 1.0),
						 0.125 * (xi + 1.0) * (eta + 1.0) * (xi + eta + 2.0 * zeta - 1.0);
			dN.col(7) << -0.125 * (eta + 1.0) * (zeta + 1.0) * (-2.0 * xi + eta + zeta - 1.0),
						 -0.125 * (xi - 1.0) * (zeta + 1.0) * (-xi + 2.0 * eta + zeta - 1.0),
						 -0.125 * (xi - 1.0) * (eta + 1.0) * (-xi + eta + 2.0 * zeta - 1.0);
			dN.col(8) << -0.5 * xi * (eta - 1.0) * (zeta - 1.0),
						 -0.25 * (xi2 - 1.0) * (zeta - 1.0),
						 -0.25 * (xi2 - 1.0) * (eta - 1.0);
			dN.col(9) << 0.25 * (eta2 - 1.0) * (zeta - 1.0),
						 0.5 * (xi + 1.0) * eta * (zeta - 1.0),
						 0.25 * (xi + 1.0) * (eta2 - 1.0);
			dN.col(10) << 0.5 * xi * (eta + 1.0) * (zeta - 1.0),
						  0.25 * (xi2 - 1.0) * (zeta - 1.0),
						  0.25 * (xi2 - 1.0) * (eta + 1.0);
			dN.col(11) << -0.25 * (eta2 - 1.0) * (zeta - 1.0),
						  -0.5 * (xi - 1.0) * eta * (zeta - 1.0),
						 - 0.25 * (xi - 1.0) * (eta2 - 1.0);
			dN.col(12) << 0.5 * xi * (eta - 1.0) * (zeta + 1.0),
						  0.25 * (xi2 - 1.0) * (zeta + 1.0),
						  0.25 * (xi2 - 1.0) * (eta - 1.0);
			dN.col(13) << -0.25 * (eta2 - 1.0) * (zeta + 1.0),
						  -0.5 * (xi + 1.0) * eta * (zeta + 1.0),
						  -0.25 * (xi + 1.0) * (eta2 - 1.0);
			dN.col(14) << -0.5 * xi * (eta + 1.0) * (zeta + 1.0),
						  -0.25 * (xi2 - 1.0) * (zeta + 1.0),
						  -0.25 * (xi2 - 1.0) * (eta + 1.0);
			dN.col(15) << 0.25 * (eta2 - 1.0) * (zeta + 1.0),
						  0.5 * (xi - 1.0) * eta * (zeta + 1.0),
						  0.25 * (xi - 1.0) * (eta2 - 1.0);
			dN.col(16) << -0.25 * (eta - 1.0) * (zeta2 - 1.0),
						  -0.25 * (xi - 1.0) * (zeta2 - 1.0),
						  -0.5 * (xi - 1.0) * (eta - 1.0) * zeta;
			dN.col(17) << 0.25 * (eta - 1.0) * (zeta2 - 1.0),
						  0.25 * (xi + 1.0) * (zeta2 - 1.0),
						  0.5 * (xi + 1.0) * (eta - 1.0) * zeta;
			dN.col(18) << -0.25 * (eta + 1.0) * (zeta2 - 1.0),
						  -0.25 * (xi + 1.0) * (zeta2 - 1.0),
						  -0.5 * (xi + 1.0) * (eta + 1.0) * zeta;
			dN.col(19) << 0.25 * (eta + 1.0) * (zeta2 - 1.0),
						  0.25 * (xi - 1.0) * (zeta2 - 1.0),
						  0.5 * (xi - 1.0) * (eta + 1.0) * zeta;

			return dN;
		}

		Hex20::QuadratureRule Hex20::quadratureRule() noexcept {
			const double a = std::sqrt(3.0 / 5.0);

			Hex20::QuadratureRule rule;

			rule.points = {
				Point(-a, -a, -a),
				Point(-a, 0.0, -a),
				Point(-a, a, -a),
				Point(0.0, -a, -a),
				Point(0.0, 0.0, -a),
				Point(0.0, a, -a),
				Point(a, -a, -a),
				Point(a, 0.0, -a),
				Point(a, a, -a),
				Point(-a, -a, 0.0),
				Point(-a, 0.0, 0.0),
				Point(-a, a, 0.0),
				Point(0.0, -a, 0.0),
				Point(0.0, 0.0, 0.0),
				Point(0.0, a, 0.0),
				Point(a, -a, 0.0),
				Point(a, 0.0, 0.0),
				Point(a, a, 0.0),
				Point(-a, -a, a),
				Point(-a, 0.0, a),
				Point(-a, a, a),
				Point(0.0, -a, a),
				Point(0.0, 0.0, a),
				Point(0.0, a, a),
				Point(a, -a, a),
				Point(a, 0.0, a),
				Point(a, a, a)
			};

			rule.weights = {
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (5/9) * (8/9) * (5/9)
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (8/9) * (5/9) * (5/9)
				320.0 / 729.0, // (8/9) * (8/9) * (5/9)
				200.0 / 729.0, // (8/9) * (5/9) * (5/9)
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (5/9) * (8/9) * (5/9)
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (5/9) * (5/9) * (8/9)
				320.0 / 729.0, // (5/9) * (8/9) * (8/9)
				200.0 / 729.0, // (5/9) * (5/9) * (8/9)
				320.0 / 729.0, // (8/9) * (5/9) * (8/9)
				512.0 / 729.0, // (8/9) * (8/9) * (8/9)
				320.0 / 729.0, // (8/9) * (5/9) * (8/9)
				200.0 / 729.0, // (5/9) * (5/9) * (8/9)
				320.0 / 729.0, // (5/9) * (8/9) * (8/9)
				200.0 / 729.0, // (5/9) * (5/9) * (8/9)
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (5/9) * (8/9) * (5/9)
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (8/9) * (5/9) * (5/9)
				320.0 / 729.0, // (8/9) * (8/9) * (5/9)
				200.0 / 729.0, // (8/9) * (5/9) * (5/9)
				125.0 / 729.0, // (5/9) * (5/9) * (5/9)
				200.0 / 729.0, // (5/9) * (8/9) * (5/9)
				125.0 / 729.0  // (5/9) * (5/9) * (5/9)
			};

			return rule;
		}

		int Hex20::gmshElementType() noexcept {
			return 17;
		}

		Hex20::NodeIndicesList Hex20::gmshNodeOrdering() noexcept {
			return {0, 1, 5, 4, 3, 2, 6, 7, 8, 16, 11, 17, 9, 12, 13, 15, 10, 19, 18, 14};
		}

	} // namespace fem

} // namespace monad
