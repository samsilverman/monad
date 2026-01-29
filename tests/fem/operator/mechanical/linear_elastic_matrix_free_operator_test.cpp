#include <tuple>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/quad4_grid.hpp"
#include "monad/grid/quad8_grid.hpp"
#include "monad/grid/hex8_grid.hpp"
#include "monad/grid/hex20_grid.hpp"
#include "monad/fem/kernel/mechanical/linear_elastic_kernel.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"
#include "monad/fem/operator/mechanical/linear_elastic_matrix_free_operator_traits.hpp"

using namespace monad;
using namespace monad::fem;
using namespace monad::fem::mechanical;

using Types = std::tuple<Quad4Grid, Quad8Grid, Hex8Grid, Hex20Grid>;

TEMPLATE_LIST_TEST_CASE("monad::fem::mechanical::LinearElasticMatrixFreeOperator: Test isSymmetric/isPSD", "[monad]", Types) {
    using Resolution = typename TestType::Resolution;
    using Size = typename TestType::Size;
    using Element = typename TestType::Element;
    using Kernel = LinearElasticKernel<Element>;
    using Operator = MatrixFreeOperator<TestType, Element, LinearElasticMatrixFreeOperatorTraits<Element::Dim>>;
    using Material = typename Kernel::Material;
    using StiffnessTensor = typename Material::MaterialTensor;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    TestType grid(resolution, size);
    grid.setDensitiesRandom(1234);

    StiffnessTensor C = StiffnessTensor::Random();
    // Make PSD
    C = C.transpose() * C;
    // Make PD
    C += StiffnessTensor::Identity();

    const Material material(C);

    const auto nodes = grid.elementNodes(0);
    const auto elementKReference = Kernel::lhs(material, nodes);

    const Operator K(grid, elementKReference);

    REQUIRE(K.isSymmetric());
    REQUIRE(K.isPSD());
}
