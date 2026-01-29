#include <stdexcept>
#include <string>
#include "monad/grid/hex8_grid.hpp"

namespace monad {

    std::size_t Hex8Grid::numNodes() const noexcept {
        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        return (nx + 1) * (ny + 1) * (nz + 1);
    }

    std::size_t Hex8Grid::numPeriodicNodes() const noexcept {
        return numElements();
    }

    Hex8Grid::Point Hex8Grid::node(std::size_t index) const {
        if (index >= numNodes()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numNodes()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        const std::size_t nodesPerPlane = (nx + 1) * (ny + 1);
        const std::size_t indexInPlane = index % nodesPerPlane;

        const double lx = size_[0];
        const double ly = size_[1];
        const double lz = size_[2];

        const double dx = lx / static_cast<double>(nx);
        const double dy = ly / static_cast<double>(ny);
        const double dz = lz / static_cast<double>(nz);

        const std::size_t i = indexInPlane % (nx + 1);
        const std::size_t j = indexInPlane / (nx + 1);
        const std::size_t k = index / nodesPerPlane;

        const double x = static_cast<double>(i) * dx;
        const double y = static_cast<double>(j) * dy;
        const double z = static_cast<double>(k) * dz;

        return Point(x, y, z);
    }

    Hex8Grid::ElementList Hex8Grid::element(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const std::size_t elementsPerPlane = nx * ny;

        const std::size_t i = index % nx;
        const std::size_t j = (index / nx) % ny;
        const std::size_t k = index / elementsPerPlane;

        auto nodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            return k * ((nx + 1) * (ny + 1)) + j * (nx + 1) + i;
        };

        // Bottom face
        const std::size_t n1 = nodeIndex(i, j, k);
        const std::size_t n2 = nodeIndex(i + 1, j, k);
        const std::size_t n3 = nodeIndex(i + 1, j + 1, k);
        const std::size_t n4 = nodeIndex(i, j + 1, k);

        // Top face
        const std::size_t n5 = nodeIndex(i, j, k + 1);
        const std::size_t n6 = nodeIndex(i + 1, j, k + 1);
        const std::size_t n7 = nodeIndex(i + 1, j + 1, k + 1);
        const std::size_t n8 = nodeIndex(i, j + 1, k + 1);

        return {n1, n2, n3, n4, n5, n6, n7, n8};
    }

    Hex8Grid::ElementList Hex8Grid::periodicElement(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];
        const std::size_t nz = resolution_[2];

        const std::size_t elementsPerPlane = nx * ny;

        const std::size_t i = index % nx;
        const std::size_t j = (index / nx) % ny;
        const std::size_t k = index / elementsPerPlane;

        auto nodeIndex = [&](std::size_t i, std::size_t j, std::size_t k) {
            const std::size_t iP = i % nx;
            const std::size_t jP = j % ny;
            const std::size_t kP = k % nz;

            return kP * (nx * ny) + jP * nx + iP;
        };

        // Bottom face
        const std::size_t n1 = nodeIndex(i, j, k);
        const std::size_t n2 = nodeIndex(i + 1, j, k);
        const std::size_t n3 = nodeIndex(i + 1, j + 1, k);
        const std::size_t n4 = nodeIndex(i, j + 1, k);

        // Top face
        const std::size_t n5 = nodeIndex(i, j, k + 1);
        const std::size_t n6 = nodeIndex(i + 1, j, k + 1);
        const std::size_t n7 = nodeIndex(i + 1, j + 1, k + 1);
        const std::size_t n8 = nodeIndex(i, j + 1, k + 1);

        return {n1, n2, n3, n4, n5, n6, n7, n8};
    }

} // namespace monad
