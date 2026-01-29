#include <stdexcept>
#include <string>
#include "monad/grid/hex20_grid.hpp"

namespace monad {

    std::size_t Hex20Grid::numNodes() const noexcept {
        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        const std::size_t numCornerNodes = (nx + 1) * (ny + 1) * (nz + 1);
        const std::size_t numXMidNodes = nx * (ny + 1) * (nz + 1);
        const std::size_t numYMidNodes = (nx + 1) * ny * (nz + 1);
        const std::size_t numZMidNodes = (nx + 1) * (ny + 1) * nz;

        return numCornerNodes + numXMidNodes + numYMidNodes + numZMidNodes;
    }

    std::size_t Hex20Grid::numPeriodicNodes() const noexcept {
        return 4 * numElements();
    }

    Hex20Grid::Point Hex20Grid::node(std::size_t index) const {
        if (index >= numNodes()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numNodes()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        const std::size_t numCornerNodes = (nx + 1) * (ny + 1) * (nz + 1);
        const std::size_t numXMidNodes = nx * (ny + 1) * (nz + 1);
        const std::size_t numYMidNodes = (nx + 1) * ny * (nz + 1);

        const double lx = size_[0];
        const double ly = size_[1];
        const double lz = size_[2];

        const double dx = lx / static_cast<double>(nx);
        const double dy = ly / static_cast<double>(ny);
        const double dz = lz / static_cast<double>(nz);

        double x;
        double y;
        double z;

        if (index < numCornerNodes) {
            const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
            const std::size_t indexInPlane = index % nodesPerPlane;

            const std::size_t i = indexInPlane % (nx + 1);
            const std::size_t j = indexInPlane / (nx + 1);
            const std::size_t k = index / nodesPerPlane;

            x = static_cast<double>(i) * dx;
            y = static_cast<double>(j) * dy;
            z = static_cast<double>(k) * dz;
        }
        else if (index < numCornerNodes + numXMidNodes) {
            // x-midpoints: grid size (i: nx, j: ny+1, k: nz+1)
            index -= numCornerNodes;

            std::size_t nodesPerPlane = nx * (ny + 1);
            std::size_t indexInPlane = index % nodesPerPlane;

            std::size_t i = indexInPlane % nx;
            std::size_t j = indexInPlane / nx;
            std::size_t k = index / nodesPerPlane;

            x = (static_cast<double>(i) + 0.5) * dx;
            y = static_cast<double>(j) * dy;
            z = static_cast<double>(k) * dz;
        }
        else if (index < numCornerNodes + numXMidNodes + numYMidNodes) {
            // y-midpoints: grid size (i: nx+1, j: ny, k: nz+1)
            index -= numCornerNodes + numXMidNodes;

            const std::size_t nodesPerPlane = (nx + 1) * ny;
            const std::size_t indexInPlane = index % nodesPerPlane;

            const std::size_t i = indexInPlane % (nx + 1);
            const std::size_t j = indexInPlane / (nx + 1);
            const std::size_t k = index / nodesPerPlane;

            x = static_cast<double>(i) * dx;
            y = (static_cast<double>(j) + 0.5) * dy;
            z = static_cast<double>(k) * dz;
        }
        else {
            // z-midpoints: grid size (i: nx+1, j: ny+1, k: nz)
            index -= numCornerNodes + numXMidNodes + numYMidNodes;

            const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
            const std::size_t indexInPlane = index % nodesPerPlane;

            const std::size_t i = indexInPlane % (nx + 1);
            const std::size_t j = indexInPlane / (nx + 1);
            const std::size_t k = index / nodesPerPlane;

            x = static_cast<double>(i) * dx;
            y = static_cast<double>(j) * dy;
            z = (static_cast<double>(k) + 0.5) * dz;
        }

        return Point(x, y, z);
    }

    Hex20Grid::ElementList Hex20Grid::element(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        const std::size_t elementsPerPlane = nx * ny;
        const std::size_t numCornerNodes = (nx + 1) * (ny + 1) * (nz + 1);
        const std::size_t numXMidNodes = nx * (ny + 1) * (nz + 1);
        const std::size_t numYMidNodes = (nx + 1) * ny * (nz + 1);

        const std::size_t i = index % nx;
        const std::size_t j = (index / nx) % ny;
        const std::size_t k = index / elementsPerPlane;

        auto cornerNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            return k * ((nx + 1) * (ny + 1)) + j * (nx + 1) + i;
        };

        auto xMidNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            return numCornerNodes + k * nx * (ny + 1) + j * nx + i;
        };

        auto yMidNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            return numCornerNodes + numXMidNodes + k * (nx + 1) * ny + j * (nx + 1) + i;
        };

        auto zMidNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            return numCornerNodes + numXMidNodes + numYMidNodes + k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
        };

        // Bottom face corners
        const std::size_t n1 = cornerNodeIndex(i, j, k);
        const std::size_t n2 = cornerNodeIndex(i + 1, j, k);
        const std::size_t n3 = cornerNodeIndex(i + 1, j + 1, k);
        const std::size_t n4 = cornerNodeIndex(i, j + 1, k);

        // Top face corners
        const std::size_t n5 = cornerNodeIndex(i, j, k + 1);
        const std::size_t n6 = cornerNodeIndex(i + 1, j, k + 1);
        const std::size_t n7 = cornerNodeIndex(i + 1, j + 1, k + 1);
        const std::size_t n8 = cornerNodeIndex(i, j + 1, k + 1);

        // Bottom face midpoints
        const std::size_t n9 = xMidNodeIndex(i, j, k); // between n1-n2
        const std::size_t n10 = yMidNodeIndex(i + 1, j, k); // between n2-n3
        const std::size_t n11 = xMidNodeIndex(i, j + 1, k); // between n3-n4
        const std::size_t n12 = yMidNodeIndex(i, j, k); // between n1-n4

        // Top face edge midpoints
        const std::size_t n13 = xMidNodeIndex(i, j, k + 1); // between n5-n6
        const std::size_t n14 = yMidNodeIndex(i + 1, j, k + 1); // between n6-n7
        const std::size_t n15 = xMidNodeIndex(i, j + 1, k + 1); // between n7-n8
        const std::size_t n16 = yMidNodeIndex(i, j, k + 1); // between n5-n8

        // Vertical midpoints
        const std::size_t n17 = zMidNodeIndex(i, j, k); // between n1-n5
        const std::size_t n18 = zMidNodeIndex(i + 1, j, k); // between n2-n6
        const std::size_t n19 = zMidNodeIndex(i + 1, j + 1, k); // between n3-n7
        const std::size_t n20 = zMidNodeIndex(i, j + 1, k); // between n4-n8

        return {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20};
    }

    Hex20Grid::ElementList Hex20Grid::periodicElement(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        const std::size_t elementsPerPlane = nx * ny;
        const std::size_t numElements = nx * ny * nz;

        const std::size_t i = index % nx;
        const std::size_t j = (index / nx) % ny;
        const std::size_t k = index / elementsPerPlane;

        auto cornerNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            // Periodic wrapping
            const std::size_t iP = (i + nx) % nx;
            const std::size_t jP = (j + ny) % ny;
            const std::size_t kP = (k + nz) % nz;

            return kP * nx * ny + jP * nx + iP;
        };

        auto xMidNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            // Periodic wrapping
            const std::size_t iP = (i + nx) % nx;
            const std::size_t jP = (j + ny) % ny;
            const std::size_t kP = (k + nz) % nz;

            return numElements + kP * nx * ny + jP * nx + iP;
        };

        auto yMidNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            // Periodic wrapping
            const std::size_t iP = (i + nx) % nx;
            const std::size_t jP = (j + ny) % ny;
            const std::size_t kP = (k + nz) % nz;

            return numElements + numElements + kP * nx * ny + jP * nx + iP;
        };

        auto zMidNodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            // Periodic wrapping
            const std::size_t iP = (i + nx) % nx;
            const std::size_t jP = (j + ny) % ny;
            const std::size_t kP = (k + nz) % nz;

            return numElements + numElements + numElements + kP * nx * ny + jP * nx + iP;
        };

        // Bottom face corners
        const std::size_t n1 = cornerNodeIndex(i, j, k);
        const std::size_t n2 = cornerNodeIndex(i + 1, j, k);
        const std::size_t n3 = cornerNodeIndex(i + 1, j + 1, k);
        const std::size_t n4 = cornerNodeIndex(i, j + 1, k);

        // Top face corners
        const std::size_t n5 = cornerNodeIndex(i, j, k + 1);
        const std::size_t n6 = cornerNodeIndex(i + 1, j, k + 1);
        const std::size_t n7 = cornerNodeIndex(i + 1, j + 1, k + 1);
        const std::size_t n8 = cornerNodeIndex(i, j + 1, k + 1);

        // Bottom face midpoints
        const std::size_t n9 = xMidNodeIndex(i, j, k); // between n1-n2
        const std::size_t n10 = yMidNodeIndex(i + 1, j, k); // between n2-n3
        const std::size_t n11 = xMidNodeIndex(i, j + 1, k); // between n3-n4
        const std::size_t n12 = yMidNodeIndex(i, j, k); // between n1-n4

        // Top face edge midpoints
        const std::size_t n13 = xMidNodeIndex(i, j, k + 1); // between n5-n6
        const std::size_t n14 = yMidNodeIndex(i + 1, j, k + 1); // between n6-n7
        const std::size_t n15 = xMidNodeIndex(i, j + 1, k + 1); // between n7-n8
        const std::size_t n16 = yMidNodeIndex(i, j, k + 1); // between n5-n8

        // Vertical midpoints
        const std::size_t n17 = zMidNodeIndex(i, j, k); // between n1-n5
        const std::size_t n18 = zMidNodeIndex(i + 1, j, k); // between n2-n6
        const std::size_t n19 = zMidNodeIndex(i + 1, j + 1, k); // between n3-n7
        const std::size_t n20 = zMidNodeIndex(i, j + 1, k); // between n4-n8

        return {n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20};
    }

} // namespace monad
