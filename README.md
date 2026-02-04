# Monad

[![Report a Bug](https://img.shields.io/static/v1.svg?label=🐛&message=Report%20a%20Bug&color=red)](https://github.com/samsilverman/monad/issues)
[![Request a Feature](https://img.shields.io/static/v1.svg?label=💡&message=Request%20a%20Feature&color=yellow)](https://github.com/samsilverman/monad/issues)

[![MacOS Build Status](https://github.com/samsilverman/monad/actions/workflows/macos-build.yml/badge.svg)](https://github.com/samsilverman/monad/actions/workflows/macos-build.yml)
[![Ubuntu Build Status](https://github.com/samsilverman/monad/actions/workflows/ubuntu-build.yml/badge.svg)](https://github.com/samsilverman/monad/actions/workflows/ubuntu-build.yml)
[![Windows Build Status](https://github.com/samsilverman/monad/actions/workflows/windows-build.yml/badge.svg)](https://github.com/samsilverman/monad/actions/workflows/windows-build.yml)

![Teaser](https://github.com/samsilverman/monad/blob/main/assets/images/teaser.svg)

> “Each simple substance has relations that express all the others.”  
> — Gottfried Wilhelm Leibniz, *Monadology* (1714)

`monad` is a C++ library for homogenization of structured 2D and 3D grids.
It provides solvers and material models for computing effective (homogenized) material tensors from heterogeneous microstructures.

## Dependencies

- **C++17 or later**
- **Eigen 3.3 or later**
- **Catch2 v3** (optional, for testing)
- **OpenMP** (optional)

## Why Monad?

Most finite element libraries are designed for general-purpose analysis, which makes periodic unit-cell homogenization difficult to implement cleanly.
Monad is built explicitly for structured grids and periodic microstructure simulations, enabling a simple and lightweight codebase.

Monad is built on standard C++, Eigen, and optional OpenMP, avoiding the heavy external dependencies of modern high-performance computing frameworks.
While GPU-based acceleration is possible in the future, Monad is designed to run easily, prioritizing simplicity, portability, and reproducibility over squeezing out every last drop of throughput.

## Available Material Models

### Mechanical Materials

| Physics | Constitutive Law | Governing PDE | Homogenized Tensor |
| - | - | - | - |
| Linear elasticity (Hooke's law) | $\boldsymbol\sigma=\mathbf{C}:\boldsymbol\varepsilon$ | $\nabla\cdot\boldsymbol\sigma=\mathbf{0}$ | Stiffness tensor $\overline{\mathbf{C}}$ |

### Transport Phenomena

| Physics | Constitutive Law | Homogenized Tensor |
| - | - | - |
| Linear dielectric (Gauss's law) | $\mathbf{D}=\boldsymbol\varepsilon\mathbf{E}$ | Permittivity tensor $\overline{\boldsymbol\varepsilon}$ |
| Linear electrical conductivity (Ohm's law) | $\mathbf{J}=\boldsymbol\sigma\mathbf{E}$ | Conductivity tensor $\overline{\boldsymbol\sigma}$ |
| Linear magnetism (current-free) | $\mathbf{B}=\boldsymbol\mu\mathbf{H}$ | Permeability tensor $\overline{\boldsymbol\mu}$ |
| Linear mass diffusion (Fick's law) | $\mathbf{J}=-\mathbf{D}\nabla c$ | Diffusivity tensor $\overline{\mathbf{D}}$ |
| Linear porosity (Darcy's law) | $\mathbf{q}=-\mathbf{K}\nabla p$ | Permeability tensor $\overline{\mathbf{K}}$ |
| Linear thermal conductivity (Fourier's law) | $\mathbf{q}=-\boldsymbol\kappa\nabla T$ | Conductivity tensor $\overline{\boldsymbol\kappa}$ |

> [!NOTE]
> In `monad`, all transport phenomena are aliases of a single underlying linear transport model, which follows the general constitutive law
>
> $$\mathbf{J}=-\mathbf{K}\nabla\phi$$
>
> and the scalar diffusion PDE
>
> $$\nabla\cdot\mathbf{J}=0.$$
>
> Only the physical interpretation of $\phi$, $\mathbf{J}$, and $\mathbf{K}$ changes.

### Multi-Physics

| Physics | Constitutive Law | Governing PDE | Homogenized Tensor |
| - | - | - | - |
| Linear piezoelectricity (stress–charge form) | $\mathbf{S}=\mathbf{c}:\mathbf{T}-\mathbf{d}^\intercal\mathbf{E}$<br>$-\mathbf{D}=-\mathbf{d}\mathbf{T}-\boldsymbol\varepsilon\mathbf{E}$ | $\nabla\cdot\mathbf{S}=\mathbf{0}$<br>$\nabla\cdot\mathbf{D}=0$ | Stiffness tensor $\overline{\mathbf{c}}$, permittivity tensor $\overline{\boldsymbol\varepsilon}$, and piezoelectric coupling tensor $\overline{\mathbf{d}}$ |

## Building

```bash
git clone https://github.com/samsilverman/monad.git
cd monad
mkdir build && cd build
cmake ..
make -j8
```

The recommended way to use `monad` is to vendor it directly (e.g., as a git submodule) and add it to your build with:

```cmake
add_subdirectory("path/to/monad")
target_link_libraries(your_target PRIVATE monad)
```

### Library Compile Flags

| CMake Flag | Default | Description |
| - | - | - |
| `MONAD_BUILD_TESTS` | `OFF` | Build test suite. |
| `MONAD_BUILD_APPS` | `OFF` | Build command-line applications. |
| `MONAD_USE_OPENMP` | `OFF` | Enable OpenMP parallelism. |

These flags are set in the cmake line:

```bash
cmake .. -DMONAD_BUILD_APPS=ON ...
```

For convenience, a `build.sh` script is included for building with compile flags:

```bash
./build.sh
```

## Minimal Example

```cpp
#include <cstddef>
#include <iostream>
#include <monad/monad.hpp>

int main() {
    // Grid resolution
    const std::size_t nx = 5;
    const std::size_t ny = 5;

    // Grid size
    const double lx = 1.0;
    const double ly = 1.0;

    monad::Quad8Grid grid({nx, ny}, {lx, ly});
    grid.setDensitiesRandom();

    // Linear elastic material
    const double E = 1.0;
    const double nu = 0.3;
    const auto planeCondition = monad::LinearElasticMaterial2d::PlaneCondition::PlaneStress;

    const monad::LinearElasticMaterial2d material(E, nu, planeCondition);

    std::cout << "Base material stiffness tensor:\n" << material.materialTensor() << std::endl;

    // Solve
    const monad::LinearElasticSolver solver(grid, material);
    const auto results = solver.solve();

    std::cout << "Homogenized stiffness tensor:\n" << results.CBar << std::endl;

    return 0;
}
```

## Documentation

Documentation for all functions is provided through Doxygen-style docstrings. [Command-line tools](https://github.com/samsilverman/monad/blob/main/apps/README.md) in `apps/` demonstrate complete example workflows.

## Applications

Build and run the provided command-line tools:

```bash
mkdir build && cd build
cmake -DMONAD_BUILD_APPS=ON ..
make -j8
./apps/<app_name>
```

See the apps [README](https://github.com/samsilverman/monad/blob/main/apps/README.md) for more details.

## Tests

Build and run the test suite:

```bash
mkdir build && cd build
cmake -DMONAD_BUILD_TESTS=ON ..
make -j8
./tests/monad_tests
```

## Contributing

Contribution guidelines are provided in [CONTRIBUTING.md](https://github.com/samsilverman/monad/blob/main/CONTRIBUTING.md).

## Roadmap

The following features are planned or under consideration for future releases:

- **Material models**: linear multi-physics models (thermoelasticity, piezomagnetism, etc.), nonlinear constitutive laws (e.g. hyperelasticity)
- **Differentiability**: gradients of homogenized tensors w.r.t. densities
- **Performance**: graph coloring, custom preconditioners, additional linear solvers

## Maintainers

- [Sam Silverman](https://github.com/samsilverman/) - [sssilver@bu.edu](mailto:sssilver@bu.edu)

## Citation

If you use Monad in your project, please consider citing the library:

```bibtex
@misc{monad,
title = {Monad},
author = {Samuel Silverman},
note = {https://github.com/samsilverman/monad},
year = {2026}
}
```

## License

Released under the [MIT License](https://github.com/samsilverman/monad/blob/main/LICENSE).
