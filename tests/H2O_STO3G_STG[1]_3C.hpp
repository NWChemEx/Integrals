#pragma once
#include "test_common_TA.hpp"

static BlockTensor corr{
        {{0, 0, 0, }, {0.4320317482681658, } },
        {{0, 0, 1, }, {0.09236922007508419, } },
        {{0, 0, 2, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.01498939895270311,
                       0.01498939895270311, } },
        {{0, 1, 0, }, {0.0923692200750842, } },
        {{0, 1, 1, }, {0.2259706497826154, } },
        {{0, 1, 2, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.0624498876847509,
                       0.0624498876847509, } },
        {{0, 2, 0, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.01498939895270311,
                       0.01498939895270311, } },
        {{0, 2, 1, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.06244988768475091,
                       0.06244988768475091, } },
        {{0, 2, 2, }, {0.2260010788313202, 0.000000000000000, 0.000000000000000, 0.02846456890651059,
                       0.02846456890651059, 0.000000000000000, 0.2260010788313201, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.2260010788313202, 0.0364327820295209, -0.0364327820295209, 0.02846456890651059,
                       0.000000000000000, 0.0364327820295209, 0.07612225412937013, 0.01914417638270967,
                       0.02846456890651059, 0.000000000000000, -0.0364327820295209, 0.01914417638270967,
                       0.07612225412937013, } },
        {{1, 0, 0, }, {1.651106781347168, } },
        {{1, 0, 1, }, {0.3814113370597547, } },
        {{1, 0, 2, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.06186067263539054,
                       0.06186067263539054, } },
        {{1, 1, 0, }, {0.3814113370597546, } },
        {{1, 1, 1, }, {1.291680373861538, } },
        {{1, 1, 2, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.4115234522687283,
                       0.4115234522687283, } },
        {{1, 2, 0, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.06186067263539053,
                       0.06186067263539053, } },
        {{1, 2, 1, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.4115234522687283,
                       0.4115234522687283, } },
        {{1, 2, 2, }, {1.291063356136088, 0.000000000000000, 0.000000000000000, 0.2035016512457414,
                       0.2035016512457414, 0.000000000000000, 1.291063356136088, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       1.291063356136088, 0.2604687717855397, -0.2604687717855397, 0.2035016512457414,
                       0.000000000000000, 0.2604687717855397, 0.6677487129306753, 0.1468365737070141,
                       0.2035016512457414, 0.000000000000000, -0.2604687717855397, 0.1468365737070141,
                       0.6677487129306753, } },
        {{2, 0, 0, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.8762474427199597,
                       0.8762474427199597, } },
        {{2, 0, 1, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.2058120840570231,
                       0.2058120840570231, } },
        {{2, 0, 2, }, {0.04403798958597818, 0.000000000000000, 0.000000000000000, 0.001458524510931188,
                       0.001458524510931188, 0.000000000000000, 0.04403798958597816, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.04403798958597818, 0.001866815751399469, -0.001866815751399469, 0.0148937824460669,
                       0.000000000000000, 0.0190630650769683, 0.03476698590254346, 0.03304865628547179,
                       0.0148937824460669, 0.000000000000000, -0.0190630650769683, 0.03304865628547179,
                       0.03476698590254346, } },
        {{2, 1, 0, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.2058120840570231,
                       0.2058120840570231, } },
        {{2, 1, 1, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.7998245551809059,
                       0.7998245551809059, } },
        {{2, 1, 2, }, {0.3452168068904369, 0.000000000000000, 0.000000000000000, 0.07061705816929115,
                       0.07061705816929115, 0.000000000000000, 0.3452168068904368, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.3452168068904369, 0.09038520471881525, -0.09038520471881525, 0.1575300389176288,
                       0.000000000000000, 0.2016281219588468, 0.4175798315807921, 0.2491307080298726,
                       0.1575300389176288, 0.000000000000000, -0.2016281219588468, 0.2491307080298726,
                       0.4175798315807921, } },
        {{2, 2, 0, }, {0.04403798958597818, 0.000000000000000, 0.000000000000000, 0.001458524510931188,
                       0.001458524510931188, 0.000000000000000, 0.04403798958597816, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.04403798958597818, 0.001866815751399469, -0.001866815751399469, 0.0148937824460669,
                       0.000000000000000, 0.0190630650769683, 0.03476698590254346, 0.03304865628547179,
                       0.0148937824460669, 0.000000000000000, -0.0190630650769683, 0.03304865628547179,
                       0.03476698590254346, } },
        {{2, 2, 1, }, {0.3452168068904369, 0.000000000000000, 0.000000000000000, 0.07061705816929115,
                       0.07061705816929115, 0.000000000000000, 0.3452168068904368, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.3452168068904369, 0.09038520471881527, -0.09038520471881527, 0.1575300389176288,
                       0.000000000000000, 0.2016281219588468, 0.4175798315807921, 0.2491307080298726,
                       0.1575300389176288, 0.000000000000000, -0.2016281219588468, 0.2491307080298726,
                       0.4175798315807921, } },
        {{2, 2, 2, }, {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.1363204531432907,
                       0.1363204531432907, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.04572648169830933, -0.04572648169830933, 0.1363204531432907,
                       0.000000000000000, 0.04572648169830933, 0.2503135844009782, 0.0531178896955531,
                       0.1363204531432907, 0.000000000000000, -0.04572648169830933, 0.0531178896955531,
                       0.2503135844009782, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.1005948094606543, 0.1005948094606543, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.1005948094606543, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.1005948094606543, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.04572648169830933, -0.04572648169830933, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, 0.000000000000000, 0.1591217081522532,
                       0.1591217081522532, 0.04572648169830933, 0.000000000000000, 0.1591217081522532,
                       0.3203849771785241, 0.000000000000000, -0.04572648169830933, 0.000000000000000,
                       0.1591217081522532, 0.000000000000000, -0.3203849771785241, 0.8029719894249459,
                       0.000000000000000, 0.03780200276131448, 0.2894539623743748, 0.170786026869631,
                       0.000000000000000, 0.7734376620290397, 0.000000000000000, 0.000000000000000,
                       0.000000000000000, 0.03780200276131448, 0.000000000000000, 0.8218217472194999,
                       0.3704820948949865, -0.0967792783441325, 0.2894539623743748, 0.000000000000000,
                       0.3704820948949865, 1.324961030886023, 0.1557568246651184, 0.170786026869631,
                       0.000000000000000, -0.0967792783441325, 0.1557568246651184, 0.4031548629489811,
                       0.8029719894249459, 0.000000000000000, -0.03780200276131448, 0.170786026869631,
                       0.2894539623743748, 0.000000000000000, 0.7734376620290397, 0.000000000000000,
                       0.000000000000000, 0.000000000000000, -0.03780200276131448, 0.000000000000000,
                       0.8218217472194999, 0.0967792783441325, -0.3704820948949865, 0.170786026869631,
                       0.000000000000000, 0.0967792783441325, 0.4031548629489811, 0.1557568246651184,
                       0.2894539623743748, 0.000000000000000, -0.3704820948949865, 0.1557568246651184,
                       1.324961030886023, } },
};