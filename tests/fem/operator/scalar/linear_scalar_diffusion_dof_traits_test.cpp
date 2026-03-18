#include <tuple>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/grid/grid_aliases.hpp"
#include "monad/field/density_field.hpp"
#include "monad/fem/kernel/scalar/linear_scalar_diffusive_kernel.hpp"
#include "monad/fem/operator/matrix_free_operator.hpp"
#include "monad/fem/operator/scalar/linear_scalar_diffusive_dof_traits.hpp"

using namespace monad;
using namespace monad::field;
using namespace monad::fem;
using namespace monad::fem::scalar;

template <class GridT, GradientConvention C>
struct TypePair {
    using Grid = GridT;
    static constexpr GradientConvention Convention = C;
};

using Types = std::tuple<
    TypePair<Quad4Grid, GradientConvention::Negative>,
    TypePair<Quad4Grid, GradientConvention::Positive>,
    TypePair<Quad8Grid, GradientConvention::Negative>,
    TypePair<Quad8Grid, GradientConvention::Positive>,
    TypePair<Hex8Grid, GradientConvention::Negative>,
    TypePair<Hex8Grid, GradientConvention::Positive>,
    TypePair<Hex20Grid, GradientConvention::Negative>,
    TypePair<Hex20Grid, GradientConvention::Positive>
>;

TEMPLATE_LIST_TEST_CASE("monad::fem::scalar: Test LinearScalarDiffusiveDofTraits", "[monad]", Types) {
    using Grid = typename TestType::Grid;
    using Resolution = typename Grid::Resolution;
    using Size = typename Grid::Size;
    using Element = typename Grid::Element;
    using Kernel = LinearScalarDiffusiveKernel<Element, TestType::Convention>;
    using Traits = LinearScalarDiffusiveDofTraits;
    using Operator = MatrixFreeOperator<Grid, Traits>;
    using Material = typename Kernel::Material;

    Resolution resolution;
    resolution.fill(2);

    Size size;
    size.fill(0.5);

    const Grid grid(resolution, size);

    DensityField<Element::Dim> densityField(resolution);
    densityField.setRandom(1234);

    const Material material(2.1);

    const auto nodes = grid.elementNodes(0);
    const auto elementKReference = Kernel::lhs(material, nodes);

    const Operator K(grid, densityField, elementKReference);

    REQUIRE(K.isSymmetric());
    REQUIRE(K.isPSD());
}
