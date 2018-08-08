set(
    INTEGRALS_SRCS nwx_libint/nwx_libint.cpp
                   LibintIntegral.cpp
)
set(
    INTEGRALS_INCLUDES LibintIntegral.hpp
                       nwx_libint/nwx_libint.hpp
)
add_nwx_library(
        NAME Integrals
        SOURCES ${INTEGRALS_SRCS}
        DEPENDS Utilities::utilities LibChemist::chemist ${TAMM_LIBRARIES}
        ${LIBINT_LIBRARY}
)
target_include_directories(Integrals PUBLIC ${TAMM_INCLUDE_DIRS}
                           ${LIBINT_INCLUDE_DIRS})
install_targets(
        TARGETS Integrals
        INCLUDES ${INTEGRALS_INCLUDES}
)
