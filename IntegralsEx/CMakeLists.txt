cmake_minimum_required(VERSION ${CMAKE_VERSION})
project(IntegralsEx-SRC VERSION 0.0.0 LANGUAGES CXX)
include(TargetMacros)
set(IntegralsEx_SRCS     TwoCTensorBuilder.cpp
                         ThreeCTensorBuilder.cpp
                         nwx_libint/nwx_libint.cpp
)
set(IntegralsEx_INCLUDES TwoCenterIntegral.hpp
                         ThreeCenterIntegral.hpp
                         FourCenterIntegral.hpp
                         TwoCTensorBuilder.hpp
                         ThreeCTensorBuilder.hpp
                         IntegralTensorBuilder.hpp
                         nwx_libint/nwx_libint.hpp
                         nwx_libint/LibInt2C.hpp
                         nwx_libint/LibInt3C.hpp
                         nwx_libint/LibInt4C.hpp
)
set(IntegralsEx_DEFINITIONS "")
set(IntegralsEx_LINKFLAGS "")
nwchemex_add_library(IntegralsEx IntegralsEx_SRCS
                                 IntegralsEx_INCLUDES
                                 IntegralsEx_DEFINTIONS
                                 IntegralsEx_LINKFLAGS
)
