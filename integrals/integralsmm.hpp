#pragma once
#include <sde/module_manager.hpp>
#include "integrals/libint_integral.hpp"

namespace integrals {
namespace libint {

void load_modules(sde::ModuleManager& mm,
                  const std::size_t tile_size = std::size_t{180},
                  const detail_::implementation_type& impl = detail_::implementation_type::direct,
                  const double screen = 0.0);

} // end namespace libint
} // end namespace integrals
