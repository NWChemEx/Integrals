set(INTEGRALS_SRCS nwx_libint/nwx_libint.cpp
                   integralsmm.cpp
                   libint_integral.cpp
                   tamm_int_functor.cpp
)
set(INTEGRALS_INCLUDES libint_integral.hpp
                       nwx_libint/nwx_libint.hpp
                       integralsmm.hpp
                       types.hpp
)
cpp_add_library(
    integrals
    SOURCES ${INTEGRALS_SRCS}
    INCLUDES ${INTEGRALS_INCLUDES}
    DEPENDS property_types::property_types libint2
)

cpp_install(TARGETS integrals)
