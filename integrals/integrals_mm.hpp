#pragma once
#include <sde/module_manager.hpp>

namespace integrals::libint {

void load_modules(sde::ModuleManager& mm,
                  std::size_t tile_size = std::size_t{180});

} // namespace integrals::libint
