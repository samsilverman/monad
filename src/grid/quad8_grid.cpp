#include <stdexcept>
#include <string>
#include "monad/grid/quad8_grid.hpp"

namespace monad {

    std::size_t Quad8Grid::numNodes() const noexcept {
        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const std::size_t numCornerNodes = (nx + 1) * (ny + 1);
        const std::size_t numXMidNodes = nx * (ny + 1);
        const std::size_t numYMidNodes = (nx + 1) * ny;

        return numCornerNodes + numXMidNodes + numYMidNodes;
    }

    std::size_t Quad8Grid::numPeriodicNodes() const noexcept {
        return 3 * numElements();
    }

    Quad8Grid::Point Quad8Grid::node(std::size_t index) const {
        if (index >= numNodes()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numNodes()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const std::size_t numCornerNodes = (nx + 1) * (ny + 1);
        const std::size_t numXMidNodes = nx * (ny + 1);

        const double lx = size_[0];
        const double ly = size_[1];

        const double dx = lx / static_cast<double>(nx);
        const double dy = ly / static_cast<double>(ny);

        double x;
        double y;

        if (index < numCornerNodes) {
            const std::size_t i = index % (nx + 1);
            const std::size_t j = index / (nx + 1);

            x = static_cast<double>(i) * dx;
            y = static_cast<double>(j) * dy;
        }
        else if (index < numCornerNodes + numXMidNodes) {
            // x-midpoints are laid out in nx columns (i) and (ny + 1) rows (j).
            index -= numCornerNodes;

            const std::size_t i = index % nx;
            const std::size_t j = index / nx;

            x = (static_cast<double>(i) + 0.5) * dx;
            y = static_cast<double>(j) * dy;
        }
        else {
            // y-midpoints are laid out in (nx + 1) columns (i) and ny rows (j).
            index -= numCornerNodes + numXMidNodes;

            const std::size_t i = index % (nx + 1);
            const std::size_t j = index / (nx + 1);

            x = static_cast<double>(i) * dx;
            y = (static_cast<double>(j) + 0.5) * dy;
        }

        return Point(x, y);
    }

    Quad8Grid::ElementList Quad8Grid::element(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const std::size_t numCornerNodes = (nx + 1) * (ny + 1);
        const std::size_t numXMidNodes = nx * (ny + 1);

        const std::size_t i = index % nx;
        const std::size_t j = index / nx;

        auto cornerNodeIndex = [&](std::size_t i, std::size_t j) {
            return j * (nx + 1) + i;
        };

        auto xMidNodeIndex = [&](std::size_t i, std::size_t j) {
            return numCornerNodes + j * nx + i;
        };

        auto yMidNodeIndex = [&](std::size_t i, std::size_t j) {
            return numCornerNodes + numXMidNodes + j * (nx + 1) + i;
        };

        // Corner nodes
        const std::size_t n1 = cornerNodeIndex(i, j);
        const std::size_t n2 = cornerNodeIndex(i + 1, j);
        const std::size_t n3 = cornerNodeIndex(i + 1, j + 1);
        const std::size_t n4 = cornerNodeIndex(i, j + 1);

        // Midpoint nodes
        const std::size_t n5 = xMidNodeIndex(i, j); // between n1-n2
        const std::size_t n6 = yMidNodeIndex(i + 1, j); // between n2-n3
        const std::size_t n7 = xMidNodeIndex(i, j + 1); // between n3-n4
        const std::size_t n8 = yMidNodeIndex(i, j); // between n1-n4

        return {n1, n2, n3, n4, n5, n6, n7, n8};
    }

    Quad8Grid::ElementList Quad8Grid::periodicElement(std::size_t index) const {
        if (index >= numElements()) {
            throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
        }

        const std::size_t nx = resolution_[0];
        const std::size_t ny = resolution_[1];

        const std::size_t numElements = nx * ny;

        const std::size_t i = index % nx;
        const std::size_t j = index / nx;

        auto cornerNodeIndex = [&](std::size_t i, std::size_t j) {
            const std::size_t iP = i % nx;
            const std::size_t jP = j % ny;

            return jP * nx + iP;
        };

        auto xMidNodeIndex = [&](std::size_t i, std::size_t j) {
            // Periodic wrapping
            const std::size_t iP = i % nx;
            const std::size_t jP = j % ny;

            return numElements + jP * nx + iP;
        };

        auto yMidNodeIndex = [&](std::size_t i, std::size_t j) {
            // Periodic wrapping
            const std::size_t iP = i % nx;
            const std::size_t jP = j % ny;

            return numElements + numElements + jP * nx + iP;
        };

        // Corner nodes
        const std::size_t n1 = cornerNodeIndex(i, j);
        const std::size_t n2 = cornerNodeIndex(i + 1, j);
        const std::size_t n3 = cornerNodeIndex(i + 1, j + 1);
        const std::size_t n4 = cornerNodeIndex(i, j + 1);

        // Midpoint nodes
        const std::size_t n5 = xMidNodeIndex(i, j); // between n1-n2
        const std::size_t n6 = yMidNodeIndex(i + 1, j); // between n2-n3
        const std::size_t n7 = xMidNodeIndex(i, j + 1); // between n3-n4
        const std::size_t n8 = yMidNodeIndex(i, j); // between n1-n4

        return {n1, n2, n3, n4, n5, n6, n7, n8};
    }

} // namespace monad
