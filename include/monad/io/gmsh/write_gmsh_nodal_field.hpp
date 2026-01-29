#pragma once

#include <ostream>
#include <stdexcept>
#include <string>
#include <Eigen/Core>

namespace monad {

    namespace io {

        namespace gmsh {

            /**
             * @brief Writes the Gmsh file `$NodeData` section containing a nodal field.
             *
             * @tparam Derived Eigen matrix type.
             *
             * @param[in,out] os Output stream.
             * @param[in] field N×1 (scalar), N×2 or N×3 (vector) field.
             * @param[in] name Optional name tag for the node field.
             *
             * @throws std::invalid_argument: if `field` is not a scalar or vector field.
             *
             * @note Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
             */
            template <typename Derived>
            void writeGmshNodalField(std::ostream &os, const Eigen::MatrixBase<Derived> &field, const std::string &name = "") {
                if (field.cols() != 1 && field.cols() != 2 && field.cols() != 3) {
                    throw std::invalid_argument("field number of columns (" + std::to_string(field.cols()) + ") must be 1, 2, or 3.");
                }

                int numStringTags = 0;
                if (!name.empty()) {
                    numStringTags = 1;
                }
                const std::string &stringTag = name;
                const int numRealTags = 0;
                // 3 tags: timestep, data dimension, number of entries.
                const int numIntegerTags = 3;

                // Section header
                os << "$NodeData\n"
                   << numStringTags << "\n";

                if (numStringTags == 1) {
                    os << "\"" << stringTag << "\"\n";
                }

                os << numRealTags << "\n"
                   << numIntegerTags << "\n"
                   << "0\n";

                if (field.cols() == 1) {
                    os << "1\n";
                }
                else {
                    os << "3\n";
                }

                os << static_cast<int>(field.rows()) << "\n";

                // Section body
                for (int i = 0; i < field.rows(); ++i) {
                    os << (i + 1);

                    const auto row = field.row(i);
                    for (int j = 0; j < field.cols(); ++j) {
                        os << " " << row(j);
                    }

                    // pad 2D vectors to 3D
                    if (field.cols() == 2) {
                        os << " 0";
                    }

                    os << "\n";
                }

                // Section footer
                os << "$EndNodeData";
            }

        } // namespace gmsh

    } // namespace io

} // namespace monad
