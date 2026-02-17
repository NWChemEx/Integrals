#pragma once
#include <array>
#include <stdexcept>
#include <type_traits>

namespace integrals::utils {

/** @brief Code factorization for looping over shell quartets.
 *
 *  @tparam AOBasisSetType The type of the AO basis sets providing the shells.
 *                         Expected to be chemist::basis_set::AOBasisSet or
 *                         const AOBasisSet.
 *
 *  Typically looping over shell quartets is done by four nested loops, one for
 *  each of the four basis sets providing shells. This causes a lot of code
 *  duplication among files, looks messy, and leads to very deep loops. This
 *  iterator wraps the looping logic into a single class which can be used in a
 *  while loop like:
 *
 *  ```
 *  ShellQuartetIterator sqi(bra0, bra1, ket0, ket1);
 *  while(sqi) {
 *      sqi.current_quartet(); // Do something with the Shell objects
 *      sqi.shell_offsets(); // Do something with the offsets of the shells
 *      ++sqi; // Move on to the next quartet
 *  }
 *  ```
 *
 *  N.b., this iterator does not keep the AO basis sets alive, so the caller
 *  must ensure that the basis sets outlive the iterator.
 */
template<typename AOBasisSetType>
class ShellQuartetIterator {
public:
    /// CV-qualified type of the AOBasisSet objects we're iterating over
    using ao_basis_set_type = AOBasisSetType;

    /// Whether the AOBasisSet type is const-qualified
    static constexpr bool is_const_v = std::is_const_v<ao_basis_set_type>;

    /// Traits class defining the types for an ao_basis_set_type object
    using ao_traits = typename ao_basis_set_type::abs_traits;

    /// Type of reference to a possibly mutable shell
    using shell_reference =
      std::conditional_t<is_const_v, typename ao_traits::const_shell_reference,
                         typename ao_traits::shell_reference>;

    /// Type used for indexing and offsets
    using size_type = typename ao_basis_set_type::size_type;

    /// Type used to return the shell offsets
    using offset_value = std::array<size_type, 4>;

    /// Type used to return the current shell quartet
    using quartet_reference = std::array<shell_reference, 4>;

    template<typename T>
    ShellQuartetIterator(T&& bra0, T&& bra1, T&& ket0, T&& ket1) :
      m_p_ao_basis_sets_{&bra0, &bra1, &ket0, &ket1},
      m_shell_offsets_{0, 0, 0, 0},
      m_done_{false} {
        if(m_p_ao_basis_sets_[0]->n_shells() == 0 ||
           m_p_ao_basis_sets_[1]->n_shells() == 0 ||
           m_p_ao_basis_sets_[2]->n_shells() == 0 ||
           m_p_ao_basis_sets_[3]->n_shells() == 0) {
            m_done_ =
              true; // If any of the basis sets have no shells, we're done
        }
    }

    /// Returns an array containing the offsets of the current shells
    decltype(auto) shell_offsets() const {
        if(m_done_) {
            throw std::out_of_range(
              "There are no valid shell quartet offsets left.");
        }
        return m_shell_offsets_;
    }

    /// Returns an array containing references to the current Shell objects
    quartet_reference current_quartet() {
        if(m_done_) {
            throw std::out_of_range("There are no valid shell quartets left.");
        }
        quartet_reference rv;
        for(size_type i = 0; i < 4; ++i) {
            const auto offset_i = m_shell_offsets_[i];
            rv[i]               = m_p_ao_basis_sets_[i]->shell(offset_i);
        }
        return rv;
    }

    // Prefix increment, sets *this to the next shell quartet and returns *this
    ShellQuartetIterator& operator++() {
        for(size_type i = 0; i < 4; ++i) {
            // Go in lexicographical order, i.e., increment from right
            const auto inverse_i = 3 - i;
            auto offset_i        = ++m_shell_offsets_[inverse_i];
            if(offset_i < m_p_ao_basis_sets_[inverse_i]->n_shells()) {
                break; // Incrementing this offset was good
            } else {
                // Increment was bad so reset it
                m_shell_offsets_[inverse_i] = 0;
            }
        }
        if(m_shell_offsets_ == offset_value{0, 0, 0, 0}) {
            m_done_ = true; // We've wrapped around so we're done
        }
        return *this;
    }

    /// Converts *this to a boolean such that it is false when iteration is done
    operator bool() const noexcept { return !m_done_; }

private:
    // Pointer to an ao_basis_set_type object
    using ao_basis_set_pointer = ao_basis_set_type*;

    // The basis sets we're iterating over.
    std::array<ao_basis_set_pointer, 4> m_p_ao_basis_sets_;

    // m_shell_offsets_[i] is the offset for basis i's shell
    offset_value m_shell_offsets_;

    // Have we gone through all the quartets yet?
    bool m_done_;
};

} // namespace integrals::utils