#include <IntegralsEx/TwoCTensorBuilder.hpp>
#define CATCH_CONFIG_MAIN
#include "TestCommon.hpp"

std::vector<std::vector<double>> corr={
    { //( 0 | 0 )
    1.0464371073132783, },
    { //( 1 | 0 )
    3.4291997049832545, },
    { //( 1 | 1 )
    26.4352259268832768, },
    { //( 2 | 0 )
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000, },
    { //( 2 | 1 )
    0.0000000000000000, -0.0000000000000002, 0.0000000000000000, },
    { //( 2 | 2 )
      5.7847979494771771, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 5.7847979494771771, 0.0000000000000000,
      0.0000000000000000, 0.0000000000000000, 5.7847979494771771, },
    { //( 3 | 0 )
    2.6052624154025357, },
    { //( 3 | 1 )
    25.3420823544301896, },
    { //( 3 | 2 )
    4.2141100758044443, 3.2924421345358774, 0.0000000000000000, },
    { //( 3 | 3 )
    39.9325704220624615, },
    { //( 4 | 0 )
    2.6052624154025357, },
    { //( 4 | 1 )
    25.3420823544301896, },
    { //( 4 | 2 )
    -4.2141100758044443, 3.2924421345358774, 0.0000000000000000, },
    { //( 4 | 3 )
    26.6712891938732142, },
    { //( 4 | 4 )
    39.9325704220624615, },
};

TEST_CASE("Testing MetricTensorBuilder"){
    
    auto atoms=make_atoms();
    auto bs=get_basis("PRIMARY",atoms);

    std::vector<LibChemist::BasisSet> basissets({bs,bs});    

    IntegralsEx::TwoCTensorBuilder<nwx_libint::Metric> metric_build;
    auto metric_tensor = metric_build.compute(atoms,basissets);

    REQUIRE(metric_tensor[0](0,0) == Approx(corr[0][0]));
    REQUIRE(metric_tensor[0](1,0) == Approx(corr[1][0]));
    
    std::vector<LibChemist::BasisSet> badsets({bs,bs,bs});  
    REQUIRE_THROWS_AS(metric_build.compute(atoms,badsets), std::length_error);  
}
