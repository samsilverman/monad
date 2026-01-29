#include "monad/io/gmsh/write_gmsh_header.hpp"

namespace monad {

    namespace io {

        namespace gmsh {

            void writeGmshHeader(std::ostream &os) noexcept {
                os << "$MeshFormat\n"
                   << "4.1 0 8\n"
                   << "$EndMeshFormat";
            }

        } // namespace gmsh

    } // namespace io

} // namespace monad
