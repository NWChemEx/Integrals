cmake_minimum_required(VERSION ${CMAKE_VERSION})
project(Integrals-SRC VERSION 0.0.0 LANGUAGES CXX)
include(TargetMacros)
set(Integrals_SRCS     TwoCTensorBuilder.cpp
                         ThreeCTensorBuilder.cpp
                         FourCTensorBuilder.cpp
                         nwx_libint/nwx_libint.cpp
)
set(Integrals_INCLUDES TwoCenterIntegral.hpp
                         ThreeCenterIntegral.hpp
                         FourCenterIntegral.hpp
                         TwoCTensorBuilder.hpp
                         ThreeCTensorBuilder.hpp
                         FourCTensorBuilder.hpp
                         IntegralTensorBuilder.hpp
                         nwx_libint/nwx_libint.hpp
                         nwx_libint/LibInt2C.hpp
                         nwx_libint/LibInt3C.hpp
                         nwx_libint/LibInt4C.hpp
)
set(Integrals_DEFINITIONS "")
set(Integrals_LINKFLAGS "")
nwchemex_add_library(Integrals Integrals_SRCS
                                 Integrals_INCLUDES
                                 Integrals_DEFINTIONS
                                 Integrals_LINKFLAGS
)
