#pragma once

namespace monad {

    /// @brief Bitmask controlling which nodal fields are stored in solver results.
    enum class FieldSave : unsigned int {
        /// @brief Do not store any nodal fields.
        None = 0,

        /// @brief Store total nodal fields.
        Total = 1 << 0,

        /// @brief Store macroscopic nodal fields.
        Macro = 1 << 1,

        /// @brief Store microscopic nodal fields.
        Micro = 1 << 2,

        /// @brief Store all nodal fields.
        All = Total | Macro | Micro
    };

    /// @brief Bitwise OR for field-save flags.
    FieldSave operator|(FieldSave a, FieldSave b) noexcept;

    /// @brief Bitwise AND for field-save flags.
    FieldSave operator&(FieldSave a, FieldSave b) noexcept;

    /**
     * @brief Returns `true` if a field-save flag is enabled, `false` otherwise.
     *
     * @param[in] flags Combined field-save flags.
     * @param[in] bit Flag to test.
     *
     * @returns `true` if `bit` is enabled in `flags`, `false` otherwise.
     */
    bool wants(FieldSave flags, FieldSave bit) noexcept;

    /// @brief Options controlling the iterative solves and result storage.
    struct SolverOptions {
        /// @brief Maximum number of solver iterations.
        int maxIterations = 1000;

        /// @brief Relative residual tolerance.
        double tolerance = 1e-6;

        /// @brief Nodal fields to store in the results.
        FieldSave fields = FieldSave::None;

        /// @brief Equality comparison.
        bool operator==(const SolverOptions &other) const;

        /// @brief Inequality comparison.
        bool operator!=(const SolverOptions &other) const;
    };

} // namespace monad
