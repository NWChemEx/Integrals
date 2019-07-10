#pragma once
#include "integrals/libint_integral.hpp"
#include <sde/module_manager.hpp>

namespace integrals::libint {

void load_modules(sde::ModuleManager& mm,
                  std::size_t tile_size = std::size_t{180},
                  const detail_::implementation_type& impl =
                    detail_::implementation_type::direct,
                  const double screen = 0.0, const double stg_exponent = 1.0);

} // namespace integrals::libint
