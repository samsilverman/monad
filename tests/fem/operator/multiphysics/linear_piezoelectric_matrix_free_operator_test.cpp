#include <tuple>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/fem/kernel/multiphysics/linear_piezoelectric_kernel.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"
#include "monad/fem/operator/multiphysics/linear_piezoelectric_matrix_free_operator_traits.hpp"
#include "monad/material/mechanical/linear_elastic_material_2d.hpp"
#include "monad/material/mechanical/linear_elastic_material.hpp"
#include "monad/material/transport/linear_transport_material.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::multiphysics;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::fem::multiphysics::LinearPiezoelectricMatrixFreeOperator: Test isSymmetric/isPSD", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Element = typename TestType::Element;
    using Kernel = LinearPiezoelectricKernel<Element>;
    using Operator = MatrixFreeOperator<TestType, Element, LinearPiezoelectricMatrixFreeOperatorTraits<Element::Dim>>;
    using MechanicalMaterial = LinearElasticMaterial<Element::Dim>;
    using ElectricalMaterial = LinearTransportMaterial<Element::Dim>;
    using Material = typename Kernel::Material;
    using StiffnessTensor = typename MechanicalMaterial::MaterialTensor;
    using CouplingTensor = typename Material::CouplingTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);
    grid.setDensitiesRandom(1234);

    StiffnessTensor c = StiffnessTensor::Random();
    // Make PSD
    c = c.transpose() * c;
    // Make PD
    c += StiffnessTensor::Identity();

    const MechanicalMaterial elasticMaterial(c);
    const ElectricalMaterial dielectricMaterial(2.1);

    const CouplingTensor d = 0.1 * CouplingTensor::Random();

    const Material material(elasticMaterial, dielectricMaterial, d);

    const auto nodes = grid.elementNodes(0);
    const auto elementKReference = Kernel::lhs(material, nodes);

    const Operator K(grid, elementKReference);

    REQUIRE(K.isSymmetric());
    REQUIRE(!K.isPSD());
}
