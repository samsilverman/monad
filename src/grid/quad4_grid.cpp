#include <stdexcept>
#include <string>
#include "monad/grid/quad4_grid.hpp"

namespace monad {

    std::size_t Quad4Grid::numNodes() const noexcept {
        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        return (nx + 1) * (ny + 1);
    }

    std::size_t Quad4Grid::numPeriodicNodes() const noexcept {
        return numElements();
    }

    Quad4Grid::Point Quad4Grid::node(std::size_t index) const {
        if (index >= numNodes()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numNodes()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const double lx = size_[0];
        const double ly = size_[1];

        const double dx = lx / static_cast<double>(nx);
        const double dy = ly / static_cast<double>(ny);

        const std::size_t i = index % (nx + 1);
        const std::size_t j = index / (nx + 1);

        const double x = static_cast<double>(i) * dx;
        const double y = static_cast<double>(j) * dy;

        return Point(x, y);
    }

    Quad4Grid::ElementList Quad4Grid::element(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];

        const std::size_t i = index % nx;
        const std::size_t j = index / nx;

        auto nodeIndex = [&](std::size_t i, std::size_t j) {
            return j * (nx + 1) + i;
        };

        const std::size_t n1 = nodeIndex(i, j);
        const std::size_t n2 = nodeIndex(i + 1, j);
        const std::size_t n3 = nodeIndex(i + 1, j + 1);
        const std::size_t n4 = nodeIndex(i, j + 1);

        return {n1, n2, n3, n4};
    }

    Quad4Grid::ElementList Quad4Grid::periodicElement(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const std::size_t i = index % nx;
        const std::size_t j = index / nx;

        auto nodeIndex = [&](std::size_t i, std::size_t j) {
            const std::size_t iP = i % nx;
            const std::size_t jP = j % ny;

            return jP * nx + iP;
        };

        const std::size_t n1 = nodeIndex(i, j);
        const std::size_t n2 = nodeIndex(i + 1, j);
        const std::size_t n3 = nodeIndex(i + 1, j + 1);
        const std::size_t n4 = nodeIndex(i, j + 1);

        return {n1, n2, n3, n4};
    }

} // namespace monad
