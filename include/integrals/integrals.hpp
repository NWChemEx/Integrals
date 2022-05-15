/** @file integrals.hpp
 *
 *  This is a convenience header which brings in the entire public API of the
 *  integrals library. It should not be included anywhere in the integrals
 *  source files or public header files (it's okay to use it in the tests if
 *  your unit test also needs most of the headers included by it).
 */
#pragma once
#include "integrals/factory.hpp"
#include "integrals/integralsmm.hpp"
#include "simde/tensor_representation/tensor_representation.hpp"
#include "simde/types.hpp"
