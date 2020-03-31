#pragma once
#include "test_common_TA.hpp"

static BlockTensor corr{
        {{0, 0, 0, }, {  0.4320317482681658,  0.0923692200750842,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0149893989527031,  0.0149893989527031,  0.0923692200750842,
                              0.2259706497826154,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0624498876847509,  0.0624498876847509,  0.0000000000000000,  0.0000000000000000,
                              0.2260010788313202,  0.0000000000000000,  0.0000000000000000,  0.0284645689065106,
                              0.0284645689065106,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.2260010788313201,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.2260010788313202,  0.0364327820295209, -0.0364327820295209,  0.0149893989527031,
                              0.0624498876847509,  0.0284645689065106,  0.0000000000000000,  0.0364327820295209,
                              0.0761222541293701,  0.0191441763827097,  0.0149893989527031,  0.0624498876847509,
                              0.0284645689065106,  0.0000000000000000, -0.0364327820295209,  0.0191441763827097,
                              0.0761222541293701,  1.6511067813471678,  0.3814113370597547,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0618606726353905,  0.0618606726353905,
                              0.3814113370597546,  1.2916803738615377,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.4115234522687283,  0.4115234522687283,  0.0000000000000000,
                              0.0000000000000000,  1.2910633561360882,  0.0000000000000000,  0.0000000000000000,
                              0.2035016512457414,  0.2035016512457414,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  1.2910633561360878,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  1.2910633561360882,  0.2604687717855397, -0.2604687717855397,
                              0.0618606726353905,  0.4115234522687283,  0.2035016512457414,  0.0000000000000000,
                              0.2604687717855397,  0.6677487129306753,  0.1468365737070141,  0.0618606726353905,
                              0.4115234522687283,  0.2035016512457414,  0.0000000000000000, -0.2604687717855397,
                              0.1468365737070141,  0.6677487129306753,  0.0000000000000000,  0.0000000000000000,
                              0.0440379895859782,  0.0000000000000000,  0.0000000000000000,  0.0014585245109312,
                              0.0014585245109312,  0.0000000000000000,  0.0000000000000000,  0.3452168068904369,
                              0.0000000000000000,  0.0000000000000000,  0.0706170581692911,  0.0706170581692911,
                              0.0440379895859782,  0.3452168068904369,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.1363204531432907,  0.1363204531432907,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0457264816983093,
                              -0.0457264816983093,  0.0014585245109312,  0.0706170581692911,  0.1363204531432907,
                              0.0000000000000000,  0.0457264816983093,  0.2503135844009782,  0.0531178896955531,
                              0.0014585245109312,  0.0706170581692911,  0.1363204531432907,  0.0000000000000000,
                              -0.0457264816983093,  0.0531178896955531,  0.2503135844009782,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0440379895859782,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.3452168068904368,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0440379895859782,  0.3452168068904368,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.1005948094606543,  0.1005948094606543,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.1005948094606543,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.1005948094606543,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0440379895859782,  0.0018668157513995, -0.0018668157513995,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.3452168068904369,
                              0.0903852047188153, -0.0903852047188153,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0457264816983093,
                              -0.0457264816983093,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000,
                              0.0440379895859782,  0.3452168068904369,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.1591217081522532,  0.1591217081522532,  0.0018668157513995,
                              0.0903852047188153,  0.0457264816983093,  0.0000000000000000,  0.1591217081522532,
                              0.3203849771785241,  0.0000000000000000, -0.0018668157513995, -0.0903852047188153,
                              -0.0457264816983093,  0.0000000000000000,  0.1591217081522532,  0.0000000000000000,
                              -0.3203849771785241,  0.8762474427199597,  0.2058120840570231,  0.0148937824460669,
                              0.0000000000000000,  0.0190630650769683,  0.0347669859025435,  0.0330486562854718,
                              0.2058120840570231,  0.7998245551809059,  0.1575300389176288,  0.0000000000000000,
                              0.2016281219588468,  0.4175798315807921,  0.2491307080298726,  0.0148937824460669,
                              0.1575300389176288,  0.8029719894249459,  0.0000000000000000,  0.0378020027613145,
                              0.2894539623743748,  0.1707860268696310,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.7734376620290397,  0.0000000000000000,  0.0000000000000000,
                              0.0000000000000000,  0.0190630650769683,  0.2016281219588468,  0.0378020027613145,
                              0.0000000000000000,  0.8218217472194999,  0.3704820948949865, -0.0967792783441325,
                              0.0347669859025435,  0.4175798315807921,  0.2894539623743748,  0.0000000000000000,
                              0.3704820948949865,  1.3249610308860233,  0.1557568246651184,  0.0330486562854718,
                              0.2491307080298726,  0.1707860268696310,  0.0000000000000000, -0.0967792783441325,
                              0.1557568246651184,  0.4031548629489811,  0.8762474427199597,  0.2058120840570231,
                              0.0148937824460669,  0.0000000000000000, -0.0190630650769683,  0.0330486562854718,
                              0.0347669859025435,  0.2058120840570231,  0.7998245551809059,  0.1575300389176288,
                              0.0000000000000000, -0.2016281219588468,  0.2491307080298726,  0.4175798315807921,
                              0.0148937824460669,  0.1575300389176288,  0.8029719894249459,  0.0000000000000000,
                              -0.0378020027613145,  0.1707860268696310,  0.2894539623743748,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000,  0.7734376620290397,  0.0000000000000000,
                              0.0000000000000000,  0.0000000000000000, -0.0190630650769683, -0.2016281219588468,
                              -0.0378020027613145,  0.0000000000000000,  0.8218217472194999,  0.0967792783441325,
                              -0.3704820948949865,  0.0330486562854718,  0.2491307080298726,  0.1707860268696310,
                              0.0000000000000000,  0.0967792783441325,  0.4031548629489811,  0.1557568246651184,
                              0.0347669859025435,  0.4175798315807921,  0.2894539623743748,  0.0000000000000000,
                              -0.3704820948949865,  0.1557568246651184,  1.3249610308860233, } },
};
