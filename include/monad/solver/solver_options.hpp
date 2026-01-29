#pragma once

namespace monad {

    /// @brief Bitmask controlling which nodal fields are stored in solver results.
    enum class FieldSave : unsigned int {
        /// @brief Do not store any noal fields.
        None = 0,

        /// @brief Store total nodal fields X=X̄+X̃.
        Total = 1 << 0,

        /// @brief Store macroscopic nodal fields X̄.
        Macro = 1 << 1,

        /// @brief Store microscopic nodal fields X̃.
        Micro = 1 << 2,

        /// @brief Store all nodal fields.
        All = Total | Macro | Micro
    };

    /// @brief Bitwise OR operator.
    FieldSave operator|(FieldSave a, FieldSave b) noexcept;

    /// @brief Bitwise AND operator.
    FieldSave operator&(FieldSave a, FieldSave b) noexcept;

    /**
     * @brief Checks whether a specific field-save flag is enabled.
     *
     * @param[in] flags Combined field-save flags.
     * @param[in] bit Specific field-save option to test.
     *
     * @returns `true` if `bit` is enabled in `flags`, `false` otherwise.
     */
    bool wants(FieldSave flags, FieldSave bit) noexcept;

    /// @brief Solver options.
    struct SolverOptions {
        /// @brief Maximum number of solver iterations.
        int maxIterations = 1000;

        /// @brief Convergence tolerance for the solver's residual norm.
        double tolerance = 1e-6;

        /// @brief Controls which nodal fields are stored in the solver results.
        FieldSave fields = FieldSave::None;

        /// @brief Default solver options.
        static SolverOptions defaults();

        /// @brief Equality comparison.
        bool operator==(const SolverOptions &other) const;

        /// @brief Inequality comparison.
        bool operator!=(const SolverOptions &other) const;
    };

} // namespace monad
