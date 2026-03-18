#include <stdexcept>
#include <catch2/catch_template_test_macros.hpp>
#include "monad/material/transport/linear_transport_material.hpp"
#include "monad/field/density_field.hpp"
#include "monad/material/bounds.hpp"
#include "monad/detail/constants.hpp"

using namespace monad;
using namespace monad::material;
using namespace monad::field;

using Types = std::tuple<LinearTransportMaterial<2>, LinearTransportMaterial<3>>;

TEMPLATE_LIST_TEST_CASE("monad: Test voigtBound/reussBound", "[monad]", Types) {
    using Resolution = typename DensityField<TestType::Dim>::Resolution;

    const TestType material(2.1);

    Resolution resolution;
    resolution.fill(2);

    DensityField<TestType::Dim> densityField(resolution);

    SECTION("Solid material: reuss=voigt=K") {
        const auto K = material.materialTensor();

        densityField.setOnes();

        const auto voigt = voigtBound(material, densityField);
        const auto reuss = reussBound(material, densityField);

        REQUIRE(voigt.isApprox(K, NUMERICAL_ZERO));
        REQUIRE(reuss.isApprox(K, NUMERICAL_ZERO));
    }

    SECTION("Random material: reuss≤voigt") {
        densityField.setRandom(1234);
        
        const auto voigt = voigtBound(material, densityField);
        const auto reuss = reussBound(material, densityField);    

        REQUIRE(reuss.trace() <= voigt.trace());
    }
}
