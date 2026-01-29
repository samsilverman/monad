#include "monad/solver/solver_options.hpp"

namespace monad {

    FieldSave operator|(FieldSave a, FieldSave b) noexcept {
        return static_cast<FieldSave>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
    }

    FieldSave operator&(FieldSave a, FieldSave b) noexcept {
        return static_cast<FieldSave>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
    }

    bool wants(FieldSave flags, FieldSave bit) noexcept {
        return (flags & bit) != FieldSave::None;
    }

    SolverOptions SolverOptions::defaults() {
        return SolverOptions{};
    }

    bool SolverOptions::operator==(const SolverOptions &other) const {
        return maxIterations == other.maxIterations
            && tolerance == other.tolerance
            && fields == other.fields;
    }

    bool SolverOptions::operator!=(const SolverOptions &other) const {
        return !(*this == other);
    }

} // namespace monad
