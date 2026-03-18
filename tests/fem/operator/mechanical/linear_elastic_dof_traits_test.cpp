#include <tuple>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/field/density_field.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/fem/operator/mechanical/linear_elastic_dof_traits.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"

using namespace monad;
using namespace monad::field;
using namespace monad::fem;
using namespace monad::fem::mechanical;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::fem::mechanical: Test LinearElasticDofTraits", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Element = typename TestType::Element;
    using Kernel = LinearElasticKernel<Element>;
    using Traits = LinearElasticDofTraits<Element::Dim>;
    using Operator = MatrixFreeOperator<TestType, Traits>;
    using Material = typename Kernel::Material;
    using StiffnessTensor = typename Material::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    const TestType grid(resolution, size);

    DensityField<Element::Dim> densityField(resolution);
    densityField.setRandom(1234);

    StiffnessTensor C = StiffnessTensor::Random();
    // Make PSD
    C = C.transpose() * C;
    // Make PD
    C += StiffnessTensor::Identity();

    const Material material(C);

    const auto nodes = grid.elementNodes(0);
    const auto elementKReference = Kernel::lhs(material, nodes);

    const Operator K(grid, densityField, elementKReference);

    REQUIRE(K.isSymmetric());
    REQUIRE(K.isPSD());
}
