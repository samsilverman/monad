#pragma once

#include "monad/field/density_field.hpp"

namespace monad {
    
    /// @brief Per-element density field for 2d grids.
    using DensityField2d = field::DensityField<2>;

    /// @brief Per-element density field for 2d grids.
    using DensityField3d = field::DensityField<3>;

} // namespace monad
