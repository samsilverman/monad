#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace monad {

    namespace grid {

        /**
         * @brief Structured grid with topology-defined node and element numbering.
         *
         * @tparam TopologyT Concrete topology class (e.g. Quad4Topology).
         */
        template <class TopologyT>
        class StructuredGrid {
        public:
            using Topology = TopologyT;
            using Element = typename Topology::Element;

            static constexpr int Dim = Element::Dim;

            using Resolution = std::array<std::size_t, Dim>;
            using Size = std::array<double, Dim>;
            using Point = typename Element::Point;
            using NodesMatrix = typename Element::NodesMatrix;
            using ElementList = std::array<std::size_t, Element::NumNodes>;
            using ElementsList = std::vector<ElementList>;
            using NodesList = std::vector<Point>;

            /**
             * @brief Constructs a structured grid.
             *
             * @param[in] resolution How many cells in each dimension.
             * @param[in] size Physical lengths in each dimension.
             *
             * @throws std::invalid_argument if any entry in `resolution` is zero.
             * @throws std::invalid_argument if any entry in `size` is non-positive.
             *
             */
            StructuredGrid(const Resolution &resolution, const Size &size)
                : resolution_(resolution), size_(size) {
                for (std::size_t i = 0; i < Dim; ++i) {
                    if (resolution_[i] == 0) {
                        throw std::invalid_argument("Resolution in dimension " + std::to_string(i + 1) + " must be positive.");
                    }

                    if (size_[i] <= 0) {
                        throw std::invalid_argument("Size in dimension " + std::to_string(i + 1) + " (" + std::to_string(size_[i]) + ") must be positive.");
                    }
                }
            }

            /// @brief How many cells in each dimension.
            const Resolution &resolution() const noexcept {
                return resolution_;
            }

            /// @brief Physical lengths in each dimension.
            const Size &size() const noexcept {
                return size_;
            }

            /// @brief Number of elements.
            std::size_t numElements() const noexcept {
                std::size_t total = 1;

                for (std::size_t n : resolution_) {
                    total *= n;
                }

                return total;
            }

            /// @brief Number of nodes.
            std::size_t numNodes() const noexcept {
                return Topology::numNodes(resolution_);
            }

            /// @brief Number of periodic nodes.
            std::size_t numPeriodicNodes() const noexcept {
                return Topology::numPeriodicNodes(resolution_);
            }

            /**
             * @brief Coordinates for a specific node.
             *
             * @param[in] index Node index.
             *
             * @returns Coordinates for a specific node.
             *
             * @throws std::out_of_range if `index` is outside the range [0,`numNodes()`).
             */
            Point node(std::size_t index) const {
                if (index >= numNodes()) {
                    throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numNodes()) + ").");
                }

                return Topology::node(index, resolution_, size_);
            }

            /// @brief Coordinates for all nodes.
            NodesList nodes() const noexcept {
                NodesList out;
                out.reserve(numNodes());

                for (std::size_t i = 0; i < numNodes(); ++i) {
                    out.push_back(node(i));
                }

                return out;
            }

            /**
             * @brief Node indices for a specific element.
             *
             * @param[in] index Element index.
             *
             * @returns Node indices for a specific element.
             *
             * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
             */
            ElementList element(std::size_t index) const {
                if (index >= numElements()) {
                    throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
                }

                return Topology::element(index, resolution_);
            }

            /// @brief Node indices for all elements.
            ElementsList elements() const noexcept {
                ElementsList out;
                out.reserve(numElements());

                for (std::size_t i = 0; i < numElements(); ++i) {
                    out.push_back(element(i));
                }

                return out;
            }

            /**
             * @brief Periodic node indices for a specific element.
             *
             * @param[in] index Element index.
             *
             * @returns Periodic node indices for a specific element.
             *
             * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
             */
            ElementList periodicElement(std::size_t index) const {
                if (index >= numElements()) {
                    throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
                }

                return Topology::periodicElement(index, resolution_);
            }

            /// @brief Periodic node indices for all elements.
            ElementsList periodicElements() const noexcept {
                ElementsList out;
                out.reserve(numElements());

                for (std::size_t i = 0; i < numElements(); ++i) {
                    out.push_back(periodicElement(i));
                }

                return out;
            }

            /**
             * @brief Nodal coordinates for a specific element.
             *
             * @param[in] index Element index.
             *
             * @returns Nodal coordinates for a specific element.
             *
             * @throws std::out_of_range if `index` is outside the range [0,`numElements()`).
             */
            NodesMatrix elementNodes(std::size_t index) const {
                if (index >= numElements()) {
                    throw std::out_of_range("Index (" + std::to_string(index) + ") is out of range [0," + std::to_string(numElements()) + ").");
                }

                NodesMatrix out;

                int row = 0;
                for (std::size_t nodeIndex : element(index)) {
                    out.row(row++) = node(nodeIndex);
                }

                return out;
            }

            /// @brief Grid area (2D) or volume (3D).
            double measure() const noexcept {
                const auto referenceElementNodes = elementNodes(0);
                const double referenceElementMeasure = Element::measure(referenceElementNodes);

                return referenceElementMeasure * static_cast<double>(numElements());
            }

            /// @brief Grid area.
            template <int D = Dim>
            std::enable_if_t<D == 2, double> area() const noexcept {
                return measure();
            }

            /// @brief Grid volume.
            template <int D = Dim>
            std::enable_if_t<D == 3, double> volume() const noexcept {
                return measure();
            }

            /// @brief Equality comparison.
            bool operator==(const StructuredGrid<Topology> &other) const noexcept {
                return resolution_ == other.resolution_ && size_ == other.size_;
            }

            /// @brief Inequality comparison.
            bool operator!=(const StructuredGrid<Topology> &other) const noexcept {
                return !(*this == other);
            }

        private:
            /// @brief How many cells in each dimension.
            Resolution resolution_;

            /// @brief Physical lengths in each dimension.
            Size size_;
        };

    } // namespace grid

} // namespace monad
