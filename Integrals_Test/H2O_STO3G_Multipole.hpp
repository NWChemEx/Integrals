#pragma once
#include "TestCommon.hpp"

//This file meant only for inclusion in multipole tests

static BlockTensor corr{
        {{0, 0, 0, }, {1.0000000000000004,}},
        {{0, 0, 1, }, {0.2367039365108476,}},
        {{0, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{0, 0, 3, }, {0.0384055905135491,}},
        {{0, 0, 4, }, {0.0384055905135491,}},
        {{0, 1, 0, }, {0.2367039365108476,}},
        {{0, 1, 1, }, {1.0000000000000002,}},
        {{0, 1, 2, }, {-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{0, 1, 3, }, {0.3861387813310929,}},
        {{0, 1, 4, }, {0.3861387813310929,}},
        {{0, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{0, 2, 1, }, {-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{0, 2, 2, }, {1.0000000000000007,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000002,0.0000000000000000,0.0000000000000000,0.0000000000000000,1.0000000000000007,}},
        {{0, 2, 3, }, {0.2097276494226498,0.0000000000000000,0.2684376412681763,}},
        {{0, 2, 4, }, {0.2097276494226498,0.0000000000000000,-0.2684376412681763,}},
        {{0, 3, 0, }, {0.0384055905135491,}},
        {{0, 3, 1, }, {0.3861387813310928,}},
        {{0, 3, 2, }, {0.2097276494226498,0.0000000000000000,0.2684376412681763,}},
        {{0, 3, 3, }, {1.0000000000000002,}},
        {{0, 3, 4, }, {0.1817608668218930,}},
        {{0, 4, 0, }, {0.0384055905135491,}},
        {{0, 4, 1, }, {0.3861387813310928,}},
        {{0, 4, 2, }, {0.2097276494226498,0.0000000000000000,-0.2684376412681763,}},
        {{0, 4, 3, }, {0.1817608668218930,}},
        {{0, 4, 4, }, {1.0000000000000002,}},
        {{1, 0, 0, }, {0.0000000000000000,}},
        {{1, 0, 1, }, {0.0000000000000000,}},
        {{1, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0507919295879124,}},
        {{1, 0, 3, }, {0.0022296487508863,}},
        {{1, 0, 4, }, {-0.0022296487508863,}},
        {{1, 1, 0, }, {0.0000000000000000,}},
        {{1, 1, 1, }, {0.0000000000000000,}},
        {{1, 1, 2, }, {0.0000000000000000,0.0000000000000000,0.6411728443249254,}},
        {{1, 1, 3, }, {0.2627419225990110,}},
        {{1, 1, 4, }, {-0.2627419225990110,}},
        {{1, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0507919295879124,}},
        {{1, 2, 1, }, {0.0000000000000000,0.0000000000000000,0.6411728443249254,}},
        {{1, 2, 2, }, {0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{1, 2, 3, }, {0.1473996438105677,0.0000000000000000,0.4376297982033156,}},
        {{1, 2, 4, }, {-0.1473996438105677,0.0000000000000000,0.4376297982033156,}},
        {{1, 3, 0, }, {0.0022296487508863,}},
        {{1, 3, 1, }, {0.2627419225990111,}},
        {{1, 3, 2, }, {0.1473996438105677,0.0000000000000000,0.4376297982033157,}},
        {{1, 3, 3, }, {1.6380335020342405,}},
        {{1, 3, 4, }, {-0.0000000000000000,}},
        {{1, 4, 0, }, {-0.0022296487508863,}},
        {{1, 4, 1, }, {-0.2627419225990111,}},
        {{1, 4, 2, }, {-0.1473996438105677,0.0000000000000000,0.4376297982033157,}},
        {{1, 4, 3, }, {0.0000000000000000,}},
        {{1, 4, 4, }, {-1.6380335020342405,}},
        {{2, 0, 0, }, {-0.1432223429807861,}},
        {{2, 0, 1, }, {-0.0339012923798588,}},
        {{2, 0, 2, }, {0.0507919295879124,0.0000000000000000,0.0000000000000000,}},
        {{2, 0, 3, }, {-0.0037585363454903,}},
        {{2, 0, 4, }, {-0.0037585363454903,}},
        {{2, 1, 0, }, {-0.0339012923798588,}},
        {{2, 1, 1, }, {-0.1432223429807860,}},
        {{2, 1, 2, }, {0.6411728443249253,0.0000000000000000,0.0000000000000000,}},
        {{2, 1, 3, }, {0.1499739403064126,}},
        {{2, 1, 4, }, {0.1499739403064126,}},
        {{2, 2, 0, }, {0.0507919295879124,0.0000000000000000,0.0000000000000000,}},
        {{2, 2, 1, }, {0.6411728443249254,0.0000000000000000,0.0000000000000000,}},
        {{2, 2, 2, }, {-0.1432223429807861,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1432223429807860,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1432223429807861,}},
        {{2, 2, 3, }, {0.3340921027631434,0.0000000000000000,0.1089533758839038,}},
        {{2, 2, 4, }, {0.3340921027631434,0.0000000000000000,-0.1089533758839038,}},
        {{2, 3, 0, }, {-0.0037585363454903,}},
        {{2, 3, 1, }, {0.1499739403064126,}},
        {{2, 3, 2, }, {0.3340921027631434,0.0000000000000000,0.1089533758839038,}},
        {{2, 3, 3, }, {1.1365568803584103,}},
        {{2, 3, 4, }, {0.2065815637663311,}},
        {{2, 4, 0, }, {-0.0037585363454903,}},
        {{2, 4, 1, }, {0.1499739403064126,}},
        {{2, 4, 2, }, {0.3340921027631434,0.0000000000000000,-0.1089533758839038,}},
        {{2, 4, 3, }, {0.2065815637663311,}},
        {{2, 4, 4, }, {1.1365568803584103,}},
        {{3, 0, 0, }, {0.0000000000000000,}},
        {{3, 0, 1, }, {0.0000000000000000,}},
        {{3, 0, 2, }, {0.0000000000000000,0.0507919295879124,0.0000000000000000,}},
        {{3, 0, 3, }, {0.0000000000000000,}},
        {{3, 0, 4, }, {0.0000000000000000,}},
        {{3, 1, 0, }, {0.0000000000000000,}},
        {{3, 1, 1, }, {0.0000000000000000,}},
        {{3, 1, 2, }, {0.0000000000000000,0.6411728443249253,0.0000000000000000,}},
        {{3, 1, 3, }, {0.0000000000000000,}},
        {{3, 1, 4, }, {0.0000000000000000,}},
        {{3, 2, 0, }, {0.0000000000000000,0.0507919295879124,0.0000000000000000,}},
        {{3, 2, 1, }, {0.0000000000000000,0.6411728443249253,0.0000000000000000,}},
        {{3, 2, 2, }, {0.0000000000000000,-0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{3, 2, 3, }, {0.0000000000000000,0.2489679178208190,0.0000000000000000,}},
        {{3, 2, 4, }, {0.0000000000000000,0.2489679178208190,0.0000000000000000,}},
        {{3, 3, 0, }, {0.0000000000000000,}},
        {{3, 3, 1, }, {0.0000000000000000,}},
        {{3, 3, 2, }, {0.0000000000000000,0.2489679178208190,0.0000000000000000,}},
        {{3, 3, 3, }, {0.0000000000000000,}},
        {{3, 3, 4, }, {0.0000000000000000,}},
        {{3, 4, 0, }, {0.0000000000000000,}},
        {{3, 4, 1, }, {0.0000000000000000,}},
        {{3, 4, 2, }, {0.0000000000000000,0.2489679178208190,0.0000000000000000,}},
        {{3, 4, 3, }, {0.0000000000000000,}},
        {{3, 4, 4, }, {0.0000000000000000,}},
        {{4, 0, 0, }, {0.0170208494985677,}},
        {{4, 0, 1, }, {0.0138640214103877,}},
        {{4, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{4, 0, 3, }, {0.0023816767849899,}},
        {{4, 0, 4, }, {0.0023816767849899,}},
        {{4, 1, 0, }, {0.0138640214103877,}},
        {{4, 1, 1, }, {0.4916643304992799,}},
        {{4, 1, 2, }, {-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{4, 1, 3, }, {0.4609027641354013,}},
        {{4, 1, 4, }, {0.4609027641354013,}},
        {{4, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{4, 2, 1, }, {-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{4, 2, 2, }, {0.2970906676295372,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.2970906676295371,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.8912720028886116,}},
        {{4, 2, 3, }, {0.2300037785707814,0.0000000000000000,0.5566600258378568,}},
        {{4, 2, 4, }, {0.2300037785707814,0.0000000000000000,-0.5566600258378568,}},
        {{4, 3, 0, }, {0.0023816767849899,}},
        {{4, 3, 1, }, {0.4609027641354014,}},
        {{4, 3, 2, }, {0.2300037785707815,0.0000000000000000,0.5566600258378570,}},
        {{4, 3, 3, }, {3.3326779899270678,}},
        {{4, 3, 4, }, {0.2668878295504835,}},
        {{4, 4, 0, }, {0.0023816767849899,}},
        {{4, 4, 1, }, {0.4609027641354014,}},
        {{4, 4, 2, }, {0.2300037785707815,0.0000000000000000,-0.5566600258378570,}},
        {{4, 4, 3, }, {0.2668878295504835,}},
        {{4, 4, 4, }, {3.3326779899270678,}},
        {{5, 0, 0, }, {0.0000000000000000,}},
        {{5, 0, 1, }, {0.0000000000000000,}},
        {{5, 0, 2, }, {0.0000000000000000,0.0000000000000000,-0.0072745391600959,}},
        {{5, 0, 3, }, {-0.0001466458948945,}},
        {{5, 0, 4, }, {0.0001466458948945,}},
        {{5, 1, 0, }, {0.0000000000000000,}},
        {{5, 1, 1, }, {0.0000000000000000,}},
        {{5, 1, 2, }, {0.0000000000000000,0.0000000000000000,-0.0918302770198706,}},
        {{5, 1, 3, }, {0.1308140832135668,}},
        {{5, 1, 4, }, {-0.1308140832135668,}},
        {{5, 2, 0, }, {0.0000000000000000,0.0000000000000000,-0.0072745391600959,}},
        {{5, 2, 1, }, {0.0000000000000000,0.0000000000000000,-0.0918302770198706,}},
        {{5, 2, 2, }, {0.0000000000000000,0.0000000000000000,0.2970906676295372,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.2970906676295372,0.0000000000000000,0.0000000000000000,}},
        {{5, 2, 3, }, {0.2096771342787588,0.0000000000000000,0.1673254135138940,}},
        {{5, 2, 4, }, {-0.2096771342787588,0.0000000000000000,0.1673254135138940,}},
        {{5, 3, 0, }, {-0.0001466458948944,}},
        {{5, 3, 1, }, {0.1308140832135669,}},
        {{5, 3, 2, }, {0.2096771342787589,0.0000000000000000,0.1673254135138940,}},
        {{5, 3, 3, }, {1.8617182469945979,}},
        {{5, 3, 4, }, {-0.0000000000000000,}},
        {{5, 4, 0, }, {0.0001466458948944,}},
        {{5, 4, 1, }, {-0.1308140832135669,}},
        {{5, 4, 2, }, {-0.2096771342787589,0.0000000000000000,0.1673254135138940,}},
        {{5, 4, 3, }, {0.0000000000000000,}},
        {{5, 4, 4, }, {-1.8617182469945979,}},
        {{6, 0, 0, }, {0.0000000000000000,}},
        {{6, 0, 1, }, {0.0000000000000000,}},
        {{6, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{6, 0, 3, }, {0.0000000000000000,}},
        {{6, 0, 4, }, {0.0000000000000000,}},
        {{6, 1, 0, }, {0.0000000000000000,}},
        {{6, 1, 1, }, {0.0000000000000000,}},
        {{6, 1, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{6, 1, 3, }, {0.0000000000000000,}},
        {{6, 1, 4, }, {0.0000000000000000,}},
        {{6, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{6, 2, 1, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{6, 2, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.2970906676295371,0.0000000000000000,0.2970906676295371,0.0000000000000000,}},
        {{6, 2, 3, }, {0.0000000000000000,0.1311351342935433,0.0000000000000000,}},
        {{6, 2, 4, }, {0.0000000000000000,-0.1311351342935433,0.0000000000000000,}},
        {{6, 3, 0, }, {0.0000000000000000,}},
        {{6, 3, 1, }, {0.0000000000000000,}},
        {{6, 3, 2, }, {0.0000000000000000,0.1311351342935434,0.0000000000000000,}},
        {{6, 3, 3, }, {0.0000000000000000,}},
        {{6, 3, 4, }, {0.0000000000000000,}},
        {{6, 4, 0, }, {0.0000000000000000,}},
        {{6, 4, 1, }, {0.0000000000000000,}},
        {{6, 4, 2, }, {0.0000000000000000,-0.1311351342935434,0.0000000000000000,}},
        {{6, 4, 3, }, {0.0000000000000000,}},
        {{6, 4, 4, }, {0.0000000000000000,}},
        {{7, 0, 0, }, {0.0375334890274736,}},
        {{7, 0, 1, }, {0.0187194439351077,}},
        {{7, 0, 2, }, {-0.0145490783201918,0.0000000000000000,0.0000000000000000,}},
        {{7, 0, 3, }, {0.0025843787967923,}},
        {{7, 0, 4, }, {0.0025843787967923,}},
        {{7, 1, 0, }, {0.0187194439351077,}},
        {{7, 1, 1, }, {0.5121769700281860,}},
        {{7, 1, 2, }, {-0.1836605540397412,0.0000000000000000,0.0000000000000000,}},
        {{7, 1, 3, }, {0.3260288416800813,}},
        {{7, 1, 4, }, {0.3260288416800813,}},
        {{7, 2, 0, }, {-0.0145490783201918,0.0000000000000000,0.0000000000000000,}},
        {{7, 2, 1, }, {-0.1836605540397412,0.0000000000000000,0.0000000000000000,}},
        {{7, 2, 2, }, {0.9117846424175176,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.3176033071584430,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.3176033071584431,}},
        {{7, 2, 3, }, {0.2852205723600124,0.0000000000000000,0.1940725765089999,}},
        {{7, 2, 4, }, {0.2852205723600124,0.0000000000000000,-0.1940725765089999,}},
        {{7, 3, 0, }, {0.0025843787967923,}},
        {{7, 3, 1, }, {0.3260288416800813,}},
        {{7, 3, 2, }, {0.2852205723600124,0.0000000000000000,0.1940725765089999,}},
        {{7, 3, 3, }, {1.9412857784305513,}},
        {{7, 3, 4, }, {0.4115238384263898,}},
        {{7, 4, 0, }, {0.0025843787967923,}},
        {{7, 4, 1, }, {0.3260288416800813,}},
        {{7, 4, 2, }, {0.2852205723600124,0.0000000000000000,-0.1940725765089999,}},
        {{7, 4, 3, }, {0.4115238384263898,}},
        {{7, 4, 4, }, {1.9412857784305513,}},
        {{8, 0, 0, }, {0.0000000000000000,}},
        {{8, 0, 1, }, {0.0000000000000000,}},
        {{8, 0, 2, }, {0.0000000000000000,-0.0072745391600959,0.0000000000000000,}},
        {{8, 0, 3, }, {0.0000000000000000,}},
        {{8, 0, 4, }, {0.0000000000000000,}},
        {{8, 1, 0, }, {0.0000000000000000,}},
        {{8, 1, 1, }, {0.0000000000000000,}},
        {{8, 1, 2, }, {0.0000000000000000,-0.0918302770198706,0.0000000000000000,}},
        {{8, 1, 3, }, {0.0000000000000000,}},
        {{8, 1, 4, }, {0.0000000000000000,}},
        {{8, 2, 0, }, {0.0000000000000000,-0.0072745391600959,0.0000000000000000,}},
        {{8, 2, 1, }, {0.0000000000000000,-0.0918302770198706,0.0000000000000000,}},
        {{8, 2, 2, }, {0.0000000000000000,0.2970906676295371,0.0000000000000000,0.2970906676295371,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{8, 2, 3, }, {0.0000000000000000,0.0667968028392567,0.0000000000000000,}},
        {{8, 2, 4, }, {0.0000000000000000,0.0667968028392567,0.0000000000000000,}},
        {{8, 3, 0, }, {0.0000000000000000,}},
        {{8, 3, 1, }, {0.0000000000000000,}},
        {{8, 3, 2, }, {0.0000000000000000,0.0667968028392567,0.0000000000000000,}},
        {{8, 3, 3, }, {0.0000000000000000,}},
        {{8, 3, 4, }, {0.0000000000000000,}},
        {{8, 4, 0, }, {0.0000000000000000,}},
        {{8, 4, 1, }, {0.0000000000000000,}},
        {{8, 4, 2, }, {0.0000000000000000,0.0667968028392567,0.0000000000000000,}},
        {{8, 4, 3, }, {0.0000000000000000,}},
        {{8, 4, 4, }, {0.0000000000000000,}},
        {{9, 0, 0, }, {0.0170208494985677,}},
        {{9, 0, 1, }, {0.0138640214103877,}},
        {{9, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{9, 0, 3, }, {0.0021606453889894,}},
        {{9, 0, 4, }, {0.0021606453889894,}},
        {{9, 1, 0, }, {0.0138640214103877,}},
        {{9, 1, 1, }, {0.4916643304992799,}},
        {{9, 1, 2, }, {-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{9, 1, 3, }, {0.2453047234440917,}},
        {{9, 1, 4, }, {0.2453047234440917,}},
        {{9, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{9, 2, 1, }, {-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{9, 2, 2, }, {0.2970906676295372,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.8912720028886113,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.2970906676295372,}},
        {{9, 2, 3, }, {0.1024545713566023,0.0000000000000000,0.1311351342935434,}},
        {{9, 2, 4, }, {0.1024545713566023,0.0000000000000000,-0.1311351342935434,}},
        {{9, 3, 0, }, {0.0021606453889894,}},
        {{9, 3, 1, }, {0.2453047234440916,}},
        {{9, 3, 2, }, {0.1024545713566023,0.0000000000000000,0.1311351342935434,}},
        {{9, 3, 3, }, {0.6495242361405097,}},
        {{9, 3, 4, }, {0.1767321407725665,}},
        {{9, 4, 0, }, {0.0021606453889894,}},
        {{9, 4, 1, }, {0.2453047234440916,}},
        {{9, 4, 2, }, {0.1024545713566023,0.0000000000000000,-0.1311351342935434,}},
        {{9, 4, 3, }, {0.1767321407725665,}},
        {{9, 4, 4, }, {0.6495242361405097,}},
        {{10, 0, 0, }, {0.0000000000000000,}},
        {{10, 0, 1, }, {0.0000000000000000,}},
        {{10, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0079573972717598,}},
        {{10, 0, 3, }, {0.0004747037210772,}},
        {{10, 0, 4, }, {-0.0004747037210772,}},
        {{10, 1, 0, }, {0.0000000000000000,}},
        {{10, 1, 1, }, {0.0000000000000000,}},
        {{10, 1, 2, }, {0.0000000000000000,0.0000000000000000,0.7937402826910493,}},
        {{10, 1, 3, }, {0.6628075106338011,}},
        {{10, 1, 4, }, {-0.6628075106338011,}},
        {{10, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0079573972717598,}},
        {{10, 2, 1, }, {0.0000000000000000,0.0000000000000000,0.7937402826910493,}},
        {{10, 2, 2, }, {0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{10, 2, 3, }, {0.3326409283934496,0.0000000000000000,1.1408189079019369,}},
        {{10, 2, 4, }, {-0.3326409283934496,0.0000000000000000,1.1408189079019369,}},
        {{10, 3, 0, }, {0.0004747037210773,}},
        {{10, 3, 1, }, {0.6628075106338012,}},
        {{10, 3, 2, }, {0.3326409283934497,0.0000000000000000,1.1408189079019373,}},
        {{10, 3, 3, }, {7.5869231173553739,}},
        {{10, 3, 4, }, {-0.0000000000000000,}},
        {{10, 4, 0, }, {-0.0004747037210773,}},
        {{10, 4, 1, }, {-0.6628075106338012,}},
        {{10, 4, 2, }, {-0.3326409283934497,0.0000000000000000,1.1408189079019373,}},
        {{10, 4, 3, }, {0.0000000000000000,}},
        {{10, 4, 4, }, {-7.5869231173553739,}},
        {{11, 0, 0, }, {-0.0024377659447082,}},
        {{11, 0, 1, }, {-0.0019856376295315,}},
        {{11, 0, 2, }, {0.0026524657572533,0.0000000000000000,0.0000000000000000,}},
        {{11, 0, 3, }, {-0.0002025748949991,}},
        {{11, 0, 4, }, {-0.0002025748949991,}},
        {{11, 1, 0, }, {-0.0019856376295315,}},
        {{11, 1, 1, }, {-0.0704173173741864,}},
        {{11, 1, 2, }, {0.2645800942303497,0.0000000000000000,0.0000000000000000,}},
        {{11, 1, 3, }, {0.2121966928439973,}},
        {{11, 1, 4, }, {0.2121966928439973,}},
        {{11, 2, 0, }, {0.0026524657572533,0.0000000000000000,0.0000000000000000,}},
        {{11, 2, 1, }, {0.2645800942303498,0.0000000000000000,0.0000000000000000,}},
        {{11, 2, 2, }, {-0.0425500214956283,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0425500214956283,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1276500644868848,}},
        {{11, 2, 3, }, {0.3575275053016133,0.0000000000000000,0.2529147752492070,}},
        {{11, 2, 4, }, {0.3575275053016133,0.0000000000000000,-0.2529147752492070,}},
        {{11, 3, 0, }, {-0.0002025748949991,}},
        {{11, 3, 1, }, {0.2121966928439974,}},
        {{11, 3, 2, }, {0.3575275053016134,0.0000000000000000,0.2529147752492070,}},
        {{11, 3, 3, }, {3.7877780994706440,}},
        {{11, 3, 4, }, {0.3033331989595247,}},
        {{11, 4, 0, }, {-0.0002025748949991,}},
        {{11, 4, 1, }, {0.2121966928439974,}},
        {{11, 4, 2, }, {0.3575275053016134,0.0000000000000000,-0.2529147752492070,}},
        {{11, 4, 3, }, {0.3033331989595247,}},
        {{11, 4, 4, }, {3.7877780994706440,}},
        {{12, 0, 0, }, {0.0000000000000000,}},
        {{12, 0, 1, }, {0.0000000000000000,}},
        {{12, 0, 2, }, {0.0000000000000000,0.0026524657572533,0.0000000000000000,}},
        {{12, 0, 3, }, {0.0000000000000000,}},
        {{12, 0, 4, }, {0.0000000000000000,}},
        {{12, 1, 0, }, {0.0000000000000000,}},
        {{12, 1, 1, }, {0.0000000000000000,}},
        {{12, 1, 2, }, {0.0000000000000000,0.2645800942303497,0.0000000000000000,}},
        {{12, 1, 3, }, {0.0000000000000000,}},
        {{12, 1, 4, }, {0.0000000000000000,}},
        {{12, 2, 0, }, {0.0000000000000000,0.0026524657572533,0.0000000000000000,}},
        {{12, 2, 1, }, {0.0000000000000000,0.2645800942303497,0.0000000000000000,}},
        {{12, 2, 2, }, {0.0000000000000000,-0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{12, 2, 3, }, {0.0000000000000000,0.2383534426303730,0.0000000000000000,}},
        {{12, 2, 4, }, {0.0000000000000000,0.2383534426303730,0.0000000000000000,}},
        {{12, 3, 0, }, {0.0000000000000000,}},
        {{12, 3, 1, }, {0.0000000000000000,}},
        {{12, 3, 2, }, {0.0000000000000000,0.2383534426303731,0.0000000000000000,}},
        {{12, 3, 3, }, {0.0000000000000000,}},
        {{12, 3, 4, }, {0.0000000000000000,}},
        {{12, 4, 0, }, {0.0000000000000000,}},
        {{12, 4, 1, }, {0.0000000000000000,}},
        {{12, 4, 2, }, {0.0000000000000000,0.2383534426303731,0.0000000000000000,}},
        {{12, 4, 3, }, {0.0000000000000000,}},
        {{12, 4, 4, }, {0.0000000000000000,}},
        {{13, 0, 0, }, {0.0000000000000000,}},
        {{13, 0, 1, }, {0.0000000000000000,}},
        {{13, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0036943422998677,}},
        {{13, 0, 3, }, {0.0001624347556796,}},
        {{13, 0, 4, }, {-0.0001624347556796,}},
        {{13, 1, 0, }, {0.0000000000000000,}},
        {{13, 1, 1, }, {0.0000000000000000,}},
        {{13, 1, 2, }, {0.0000000000000000,0.0000000000000000,0.2777322416617103,}},
        {{13, 1, 3, }, {0.2342475452223733,}},
        {{13, 1, 4, }, {-0.2342475452223733,}},
        {{13, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0036943422998677,}},
        {{13, 2, 1, }, {0.0000000000000000,0.0000000000000000,0.2777322416617103,}},
        {{13, 2, 2, }, {0.0000000000000000,0.0000000000000000,-0.0851000429912566,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0851000429912566,0.0000000000000000,0.0000000000000000,}},
        {{13, 2, 3, }, {0.2205750811472212,0.0000000000000000,0.3335627675379246,}},
        {{13, 2, 4, }, {-0.2205750811472212,0.0000000000000000,0.3335627675379246,}},
        {{13, 3, 0, }, {0.0001624347556796,}},
        {{13, 3, 1, }, {0.2342475452223734,}},
        {{13, 3, 2, }, {0.2205750811472213,0.0000000000000000,0.3335627675379247,}},
        {{13, 3, 3, }, {3.1798911420918619,}},
        {{13, 3, 4, }, {-0.0000000000000001,}},
        {{13, 4, 0, }, {-0.0001624347556796,}},
        {{13, 4, 1, }, {-0.2342475452223734,}},
        {{13, 4, 2, }, {-0.2205750811472213,0.0000000000000000,0.3335627675379247,}},
        {{13, 4, 3, }, {0.0000000000000001,}},
        {{13, 4, 4, }, {-3.1798911420918619,}},
        {{14, 0, 0, }, {0.0000000000000000,}},
        {{14, 0, 1, }, {0.0000000000000000,}},
        {{14, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{14, 0, 3, }, {0.0000000000000000,}},
        {{14, 0, 4, }, {0.0000000000000000,}},
        {{14, 1, 0, }, {0.0000000000000000,}},
        {{14, 1, 1, }, {0.0000000000000000,}},
        {{14, 1, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{14, 1, 3, }, {0.0000000000000000,}},
        {{14, 1, 4, }, {0.0000000000000000,}},
        {{14, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{14, 2, 1, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{14, 2, 2, }, {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0425500214956283,0.0000000000000000,-0.0425500214956283,0.0000000000000000,}},
        {{14, 2, 3, }, {0.0000000000000000,0.0501898941112788,0.0000000000000000,}},
        {{14, 2, 4, }, {0.0000000000000000,-0.0501898941112788,0.0000000000000000,}},
        {{14, 3, 0, }, {0.0000000000000000,}},
        {{14, 3, 1, }, {0.0000000000000000,}},
        {{14, 3, 2, }, {0.0000000000000000,0.0501898941112788,0.0000000000000000,}},
        {{14, 3, 3, }, {0.0000000000000000,}},
        {{14, 3, 4, }, {0.0000000000000000,}},
        {{14, 4, 0, }, {0.0000000000000000,}},
        {{14, 4, 1, }, {0.0000000000000000,}},
        {{14, 4, 2, }, {0.0000000000000000,-0.0501898941112788,0.0000000000000000,}},
        {{14, 4, 3, }, {0.0000000000000000,}},
        {{14, 4, 4, }, {0.0000000000000000,}},
        {{15, 0, 0, }, {0.0000000000000000,}},
        {{15, 0, 1, }, {0.0000000000000000,}},
        {{15, 0, 2, }, {0.0000000000000000,0.0000000000000000,0.0026524657572533,}},
        {{15, 0, 3, }, {0.0001486943637824,}},
        {{15, 0, 4, }, {-0.0001486943637824,}},
        {{15, 1, 0, }, {0.0000000000000000,}},
        {{15, 1, 1, }, {0.0000000000000000,}},
        {{15, 1, 2, }, {0.0000000000000000,0.0000000000000000,0.2645800942303497,}},
        {{15, 1, 3, }, {0.1533595845181822,}},
        {{15, 1, 4, }, {-0.1533595845181822,}},
        {{15, 2, 0, }, {0.0000000000000000,0.0000000000000000,0.0026524657572533,}},
        {{15, 2, 1, }, {0.0000000000000000,0.0000000000000000,0.2645800942303497,}},
        {{15, 2, 2, }, {0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{15, 2, 3, }, {0.0689713752919001,0.0000000000000000,0.2383534426303731,}},
        {{15, 2, 4, }, {-0.0689713752919001,0.0000000000000000,0.2383534426303731,}},
        {{15, 3, 0, }, {0.0001486943637824,}},
        {{15, 3, 1, }, {0.1533595845181822,}},
        {{15, 3, 2, }, {0.0689713752919001,0.0000000000000000,0.2383534426303732,}},
        {{15, 3, 3, }, {1.0639424591813540,}},
        {{15, 3, 4, }, {-0.0000000000000000,}},
        {{15, 4, 0, }, {-0.0001486943637824,}},
        {{15, 4, 1, }, {-0.1533595845181822,}},
        {{15, 4, 2, }, {-0.0689713752919001,0.0000000000000000,0.2383534426303732,}},
        {{15, 4, 3, }, {0.0000000000000000,}},
        {{15, 4, 4, }, {-1.0639424591813540,}},
        {{16, 0, 0, }, {-0.0102511661281748,}},
        {{16, 0, 1, }, {-0.0066523178787466,}},
        {{16, 0, 2, }, {0.0110830268996031,0.0000000000000000,0.0000000000000000,}},
        {{16, 0, 3, }, {-0.0006297906791155,}},
        {{16, 0, 4, }, {-0.0006297906791155,}},
        {{16, 1, 0, }, {-0.0066523178787466,}},
        {{16, 1, 1, }, {-0.2141898204166094,}},
        {{16, 1, 2, }, {0.8331967249851309,0.0000000000000000,0.0000000000000000,}},
        {{16, 1, 3, }, {0.3056910445503893,}},
        {{16, 1, 4, }, {0.3056910445503893,}},
        {{16, 2, 0, }, {0.0110830268996031,0.0000000000000000,0.0000000000000000,}},
        {{16, 2, 1, }, {0.8331967249851308,0.0000000000000000,0.0000000000000000,}},
        {{16, 2, 2, }, {-0.3858880617547047,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1305879327809350,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1305879327809351,}},
        {{16, 2, 3, }, {0.6897787577846505,0.0000000000000000,0.1927795520312844,}},
        {{16, 2, 4, }, {0.6897787577846505,0.0000000000000000,-0.1927795520312844,}},
        {{16, 3, 0, }, {-0.0006297906791155,}},
        {{16, 3, 1, }, {0.3056910445503894,}},
        {{16, 3, 2, }, {0.6897787577846505,0.0000000000000000,0.1927795520312845,}},
        {{16, 3, 3, }, {3.6828241873072498,}},
        {{16, 3, 4, }, {0.8694525111460794,}},
        {{16, 4, 0, }, {-0.0006297906791155,}},
        {{16, 4, 1, }, {0.3056910445503894,}},
        {{16, 4, 2, }, {0.6897787577846505,0.0000000000000000,-0.1927795520312845,}},
        {{16, 4, 3, }, {0.8694525111460794,}},
        {{16, 4, 4, }, {3.6828241873072498,}},
        {{17, 0, 0, }, {0.0000000000000000,}},
        {{17, 0, 1, }, {0.0000000000000000,}},
        {{17, 0, 2, }, {0.0000000000000000,0.0036943422998677,0.0000000000000000,}},
        {{17, 0, 3, }, {0.0000000000000000,}},
        {{17, 0, 4, }, {0.0000000000000000,}},
        {{17, 1, 0, }, {0.0000000000000000,}},
        {{17, 1, 1, }, {0.0000000000000000,}},
        {{17, 1, 2, }, {0.0000000000000000,0.2777322416617102,0.0000000000000000,}},
        {{17, 1, 3, }, {0.0000000000000000,}},
        {{17, 1, 4, }, {0.0000000000000000,}},
        {{17, 2, 0, }, {0.0000000000000000,0.0036943422998677,0.0000000000000000,}},
        {{17, 2, 1, }, {0.0000000000000000,0.2777322416617102,0.0000000000000000,}},
        {{17, 2, 2, }, {0.0000000000000000,-0.0851000429912565,0.0000000000000000,-0.0851000429912565,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{17, 2, 3, }, {0.0000000000000000,0.1797206727399568,0.0000000000000000,}},
        {{17, 2, 4, }, {0.0000000000000000,0.1797206727399568,0.0000000000000000,}},
        {{17, 3, 0, }, {0.0000000000000000,}},
        {{17, 3, 1, }, {0.0000000000000000,}},
        {{17, 3, 2, }, {0.0000000000000000,0.1797206727399569,0.0000000000000000,}},
        {{17, 3, 3, }, {0.0000000000000000,}},
        {{17, 3, 4, }, {0.0000000000000000,}},
        {{17, 4, 0, }, {0.0000000000000000,}},
        {{17, 4, 1, }, {0.0000000000000000,}},
        {{17, 4, 2, }, {0.0000000000000000,0.1797206727399569,0.0000000000000000,}},
        {{17, 4, 3, }, {0.0000000000000000,}},
        {{17, 4, 4, }, {0.0000000000000000,}},
        {{18, 0, 0, }, {-0.0024377659447082,}},
        {{18, 0, 1, }, {-0.0019856376295315,}},
        {{18, 0, 2, }, {0.0026524657572533,0.0000000000000000,0.0000000000000000,}},
        {{18, 0, 3, }, {-0.0001932792728918,}},
        {{18, 0, 4, }, {-0.0001932792728918,}},
        {{18, 1, 0, }, {-0.0019856376295315,}},
        {{18, 1, 1, }, {-0.0704173173741864,}},
        {{18, 1, 2, }, {0.2645800942303497,0.0000000000000000,0.0000000000000000,}},
        {{18, 1, 3, }, {0.0846851952238480,}},
        {{18, 1, 4, }, {0.0846851952238480,}},
        {{18, 2, 0, }, {0.0026524657572533,0.0000000000000000,0.0000000000000000,}},
        {{18, 2, 1, }, {0.2645800942303497,0.0000000000000000,0.0000000000000000,}},
        {{18, 2, 2, }, {-0.0425500214956283,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1276500644868848,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0425500214956283,}},
        {{18, 2, 3, }, {0.1892874673462209,0.0000000000000000,0.0501898941112788,}},
        {{18, 2, 4, }, {0.1892874673462209,0.0000000000000000,-0.0501898941112788,}},
        {{18, 3, 0, }, {-0.0001932792728918,}},
        {{18, 3, 1, }, {0.0846851952238480,}},
        {{18, 3, 2, }, {0.1892874673462209,0.0000000000000000,0.0501898941112788,}},
        {{18, 3, 3, }, {0.7382212395450372,}},
        {{18, 3, 4, }, {0.2008661305755316,}},
        {{18, 4, 0, }, {-0.0001932792728918,}},
        {{18, 4, 1, }, {0.0846851952238480,}},
        {{18, 4, 2, }, {0.1892874673462209,0.0000000000000000,-0.0501898941112788,}},
        {{18, 4, 3, }, {0.2008661305755316,}},
        {{18, 4, 4, }, {0.7382212395450372,}},
        {{19, 0, 0, }, {0.0000000000000000,}},
        {{19, 0, 1, }, {0.0000000000000000,}},
        {{19, 0, 2, }, {0.0000000000000000,0.0079573972717598,0.0000000000000000,}},
        {{19, 0, 3, }, {0.0000000000000000,}},
        {{19, 0, 4, }, {0.0000000000000000,}},
        {{19, 1, 0, }, {0.0000000000000000,}},
        {{19, 1, 1, }, {0.0000000000000000,}},
        {{19, 1, 2, }, {0.0000000000000000,0.7937402826910491,0.0000000000000000,}},
        {{19, 1, 3, }, {0.0000000000000000,}},
        {{19, 1, 4, }, {0.0000000000000000,}},
        {{19, 2, 0, }, {0.0000000000000000,0.0079573972717598,0.0000000000000000,}},
        {{19, 2, 1, }, {0.0000000000000000,0.7937402826910491,0.0000000000000000,}},
        {{19, 2, 2, }, {0.0000000000000000,-0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,}},
        {{19, 2, 3, }, {0.0000000000000000,0.4502238123051061,0.0000000000000000,}},
        {{19, 2, 4, }, {0.0000000000000000,0.4502238123051061,0.0000000000000000,}},
        {{19, 3, 0, }, {0.0000000000000000,}},
        {{19, 3, 1, }, {0.0000000000000000,}},
        {{19, 3, 2, }, {0.0000000000000000,0.4502238123051061,0.0000000000000000,}},
        {{19, 3, 3, }, {0.0000000000000000,}},
        {{19, 3, 4, }, {0.0000000000000000,}},
        {{19, 4, 0, }, {0.0000000000000000,}},
        {{19, 4, 1, }, {0.0000000000000000,}},
        {{19, 4, 2, }, {0.0000000000000000,0.4502238123051061,0.0000000000000000,}},
        {{19, 4, 3, }, {0.0000000000000000,}},
        {{19, 4, 4, }, {0.0000000000000000,}},
};
