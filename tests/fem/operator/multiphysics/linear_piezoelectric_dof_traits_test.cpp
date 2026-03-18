#include <tuple>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/field/density_field.hpp"
#include "monad/fem/kernel/multiphysics/linear_piezoelectric_kernel.hpp"
#include "monad/fem/operator/multiphysics/linear_piezoelectric_dof_traits.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"

using namespace monad;
using namespace monad::field;
using namespace monad::fem;
using namespace monad::fem::multiphysics;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::fem::multiphysics: Test LinearPiezoelectricDofTraits", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Element = typename TestType::Element;
    using Kernel = LinearPiezoelectricKernel<Element>;
    using Traits = LinearPiezoelectricDofTraits<Element::Dim>;
    using Operator = MatrixFreeOperator<TestType, Traits>;
    using Material = typename Kernel::Material;
    using StiffnessTensor = typename Material::StiffnessTensor;
    using PermittivityTensor = typename Material::PermittivityTensor;
    using CouplingTensor = typename Material::CouplingTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    const TestType grid(resolution, size);

    DensityField<Element::Dim> densityField(resolution);
    densityField.setRandom(1234);

    StiffnessTensor c = StiffnessTensor::Random();
    // Make PSD
    c = c.transpose() * c;
    // Make PD
    c += StiffnessTensor::Identity();

    PermittivityTensor epsilon = PermittivityTensor::Random();
    // Make PSD
    epsilon = epsilon.transpose() * epsilon;
    // Make PD
    epsilon += PermittivityTensor::Identity();

    const CouplingTensor d = CouplingTensor::Constant(0.01);

    const Material material(c, epsilon, d);

    const auto nodes = grid.elementNodes(0);
    const auto elementKReference = Kernel::lhs(material, nodes);

    const Operator K(grid, densityField, elementKReference);

    REQUIRE(K.isSymmetric());
    REQUIRE(!K.isPSD());
}
