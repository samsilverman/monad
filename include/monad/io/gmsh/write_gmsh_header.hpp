#pragma once

#include <ostream>

namespace monad {
    
    namespace io {

        namespace gmsh {

            /**
             * @brief Writes the Gmsh file `$MeshFormat` section.
             *
             * @param[in,out] os Output stream.
             *
             * @note Gmsh documentation: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
             */
            void writeGmshHeader(std::ostream &os) noexcept;

        } // namespace gmsh

    } // namespace io

} // namespace monad
