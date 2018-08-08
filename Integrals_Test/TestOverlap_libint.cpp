#include <Integrals/LibintIntegral.hpp>
#include "TestCommon.hpp"

using namespace Integrals;

//Computes the overlap integrals for water in STO-3G
std::vector<std::vector<double>> corr={
    {//(  0 |  0 )
    1.0000000000000004,},
    {//(  1 |  0 )
    0.2367039365108476,},
    {//(  1 |  1 )
    1.0000000000000002,},
    {//(  2 |  0 )
    0.0000000000000000,0.0000000000000000,0.0000000000000000,},
    {//(  2 |  1 )
    0.0000000000000000,-0.0000000000000000,0.0000000000000000,},
    {//(  2 |  2 )
    1.0000000000000002,0.0000000000000000,0.0000000000000000,
    0.0000000000000000,1.0000000000000002,0.0000000000000000,
    0.0000000000000000,0.0000000000000000,1.0000000000000002,},
    {//(  3 |  0 )
    0.0384055905135490,},
    {//(  3 |  1 )
    0.3861387813310925,},
    {//(  3 |  2 )
    0.2684376412681760,0.2097276494226497,0.0000000000000000,},
    {//(  3 |  3 )
    1.0000000000000002,},
    {//(  4 |  0 )
    0.0384055905135490,},
    {//(  4 |  1 )
    0.3861387813310925,},
    {//(  4 |  2 )
    -0.2684376412681760,0.2097276494226497,0.0000000000000000,},
    {//(  4 |  3 )
    0.1817608668218927,},
    {//(  4 |  4 )
    1.0000000000000002,},
};

TEST_CASE("Testing LibIntOverlap class"){
    auto [molecule, bs] = make_molecule();
    LibIntOverlap SBuilder;
    auto S = SBuilder.run(molecule, {bs, bs});
//    for(auto i = 0; i < bs.nshells(); ++i) {
//        for (auto j = 0; j <= i ; ++j) {
//            std::vector<double> buffer(bs[i].size() * bs[j].size());
//            S.get(tamm::IndexVector{i, j}, {buffer.data(), buffer.size()});
//            std::cout<<"{//( " << i <<" | " << j << ")" << std::endl;
//            for (auto &x : buffer) std::cout << x << ",";
//            std::cout << "},"<< std::endl;
//        }
//    }
    tamm::Tensor<double>::deallocate(S);
}

