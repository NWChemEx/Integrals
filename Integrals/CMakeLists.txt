set(
    INTEGRALS_SRCS nwx_libint/nwx_libint.cpp
                   LibintIntegral.cpp
)
set(
    INTEGRALS_INCLUDES LibintIntegral.hpp
                       nwx_libint/nwx_libint.hpp
)
find_package(tamm REQUIRED)
add_library(tamm INTERFACE IMPORTED)
target_link_libraries(tamm INTERFACE ${TAMM_LIBRARIES})
target_include_directories(tamm INTERFACE ${TAMM_INCLUDE_DIRS})
cpp_add_library(
    integrals
    SOURCES ${INTEGRALS_SRCS}
    INCLUDES ${INTEGRALS_INCLUDES}
    DEPENDS libchemist::chemist tamm
)

cpp_install(TARGETS integrals)
