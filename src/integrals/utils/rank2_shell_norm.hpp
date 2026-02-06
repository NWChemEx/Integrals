#pragma once
#include <array>
#include <cmath>
#include <vector>
#include <wtf/wtf.hpp>

namespace integrals::utils {

template<typename TensorType, typename IterType>
auto rank2_shell_norm(TensorType&& buffer, std::array<IterType, 4> offsets,
                      std::array<IterType, 4> naos) {
    using float_type = double;
    // Take the shell norm over the primitive integrals in this shell pair.
    float_type shell_norm = 0.0;
    std::array<IterType, 4> ao{0, 0, 0, 0};
    std::vector<IterType> abs_ao(4, 0); // Absolute indices

    using wtf::fp::float_cast;
    for(ao[0] = 0; ao[0] < naos[0]; ++ao[0]) {
        abs_ao[0] = offsets[0] + ao[0]; // Absolute index

        for(ao[1] = 0; ao[1] < naos[1]; ++ao[1]) {
            abs_ao[1] = offsets[1] + ao[1];

            for(ao[2] = 0; ao[2] < naos[2]; ++ao[2]) {
                abs_ao[2] = offsets[2] + ao[2];

                for(ao[3] = 0; ao[3] < naos[3]; ++ao[3]) {
                    abs_ao[3] = offsets[3] + ao[3];

                    const auto I_ijkl = buffer.get_elem(abs_ao);
                    auto as_float = std::fabs(float_cast<float_type>(I_ijkl));
                    shell_norm += as_float * as_float;
                } // End loop over ao[3]
            } // End loop over ao[2]
        } // End loop over ao[1]
    } // End loop over ao[0]

    return std::sqrt(shell_norm);
}
} // namespace integrals::utils