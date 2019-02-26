set(
    INTEGRALS_SRCS nwx_libint/nwx_libint.cpp
                   LibintIntegral.cpp
                   IntegralsMM.cpp
)
set(
    INTEGRALS_INCLUDES LibintIntegral.hpp
                       nwx_libint/nwx_libint.hpp
                       IntegralsMM.hpp
)
cpp_add_library(
    integrals
    SOURCES ${INTEGRALS_SRCS}
    INCLUDES ${INTEGRALS_INCLUDES}
    DEPENDS libchemist::chemist libint
)

cpp_install(TARGETS integrals)
