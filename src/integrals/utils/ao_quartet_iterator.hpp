#pragma once
#include "shell_quartet_iterator.hpp"

namespace integrals::utils {

/** @brief Iterates over AO quartets formed from the product of four bases.
 *
 *  @tparam AOBasisSetType The type of the AO basis sets providing the shells.
 *                         Expected to be chemist::basis_set::AOBasisSet or
 *                         const AOBasisSet.
 *
 *  This iterator will iterate over all AO quartets formed from the product of
 *  four basis sets. It is meant to be used in a while loop like:
 *
 *  ```
 *  AOQuartetIterator aqi(bra0, bra1, ket0, ket1);
 *  while(aqi) {
 *     aqi.shell_offsets(); // The offsets of the current shell quartet
 *
 *     // Offset of the AO quartet from start of current shell quartet
 *     aqi.relative_ao_offsets();
 *
 *     // Offset of the AO quartet from the start of the first shell quartet
 *     aqi.absolute_ao_offsets();
 *  }
 *  ```
 *
 *  This iterator loops over all AO quartets in the current shell quartet
 *  before moving on to the next shell quartet.
 *
 *
 *  N.b., this iterator does not keep the AO basis sets alive, so the caller
 *  must ensure that the basis sets outlive the iterator.
 */
template<typename AOBasisSetType>
class AOQuartetIterator {
public:
    /// Type of the AO Basis Sets we're iterating over
    using ao_basis_set_type = AOBasisSetType;

    /// How *this internally loops over the shell quartets
    using shell_quartet_iterator_type = ShellQuartetIterator<ao_basis_set_type>;

    /// Type used to return offsets of the current AO quartet
    using offset_value = typename shell_quartet_iterator_type::offset_value;

    /// Type used for indexing and offsets
    using size_type = typename offset_value::value_type;

    template<typename T>
    AOQuartetIterator(T&& bra0, T&& bra1, T&& ket0, T&& ket1) :
      m_shell_quartet_iterator_(std::forward<T>(bra0), std::forward<T>(bra1),
                                std::forward<T>(ket0), std::forward<T>(ket1)),
      m_ao_offsets_{0, 0, 0, 0},
      m_relative_offsets_{0, 0, 0, 0} {}

    /// Offset of the current shell quartet.
    decltype(auto) shell_offsets() const {
        return m_shell_quartet_iterator_.shell_offsets();
    }

    /// The shells in the current quartet
    decltype(auto) shell_quartet() {
        return m_shell_quartet_iterator_.current_quartet();
    }

    /** @brief Offsets such that first AO quartet in THIS shell is (00|00). */
    decltype(auto) relative_ao_offsets() const {
        assert_offsets_();
        return m_relative_offsets_;
    }

    /** @brief Offsets such that (00|00) is the first AO quartet in the FIRST
     *         shell quartet.
     */
    offset_value absolute_ao_offsets() const {
        assert_offsets_();

        offset_value ao_offsets;
        for(size_type i = 0; i < 4; ++i) {
            ao_offsets[i] = m_relative_offsets_[i] + m_ao_offsets_[i];
        }
        return ao_offsets;
    }

    /// Advances the iterator to the next AO quartet.
    AOQuartetIterator& operator++();

    /// Are there still more AO quartets to iterate over?
    operator bool() const noexcept {
        return static_cast<bool>(m_shell_quartet_iterator_);
    }

private:
    void increment_relative_();
    void increment_shells_();
    void assert_offsets_() const {
        if(!m_shell_quartet_iterator_) {
            throw std::out_of_range(
              "There are no valid AO quartet offsets left.");
        }
    }

    shell_quartet_iterator_type m_shell_quartet_iterator_;
    offset_value m_ao_offsets_;
    offset_value m_relative_offsets_;
};

/** @brief Prints out the current value of the iterator. */
template<typename AOBasisSetType>
std::ostream& operator<<(std::ostream& os,
                         const AOQuartetIterator<AOBasisSetType>& aqi) {
    auto shell_offsets  = aqi.shell_offsets();
    auto abs_ao_offsets = aqi.absolute_ao_offsets();
    auto rel_ao_offsets = aqi.relative_ao_offsets();
    os << "Shell offset:";
    for(std::size_t i = 0; i < 4; ++i) { os << " " << shell_offsets[i]; }
    os << "\nRelative AO offset:";
    for(std::size_t i = 0; i < 4; ++i) os << " " << rel_ao_offsets[i];
    os << "\nAbsolute AO offset: ";
    for(std::size_t i = 0; i < 4; ++i) os << " " << abs_ao_offsets[i];
    return os;
}

// -----------------------------------------------------------------------------
// -- Out of line inline definitions
// -----------------------------------------------------------------------------

template<typename AOBasisSetType>
auto AOQuartetIterator<AOBasisSetType>::operator++() -> AOQuartetIterator& {
    increment_relative_();

    // If AO quartet is (00|00) we wrapped around (or it was a (ss|ss) quartet)
    if(m_relative_offsets_ == offset_value{0, 0, 0, 0}) { increment_shells_(); }

    return *this;
}

template<typename AOBasisSetType>
void AOQuartetIterator<AOBasisSetType>::increment_relative_() {
    auto shells = m_shell_quartet_iterator_.current_quartet();
    for(size_type i = 0; i < 4; ++i) {
        // Go in lexicographical order, i.e., increment from right
        const auto inverse_i = 3 - i;
        auto offset_i        = ++m_relative_offsets_[inverse_i];
        if(offset_i < shells[inverse_i].size()) {
            return; // Incrementing this offset was good
        } else {
            // Increment was bad so reset it
            m_relative_offsets_[inverse_i] = 0;
        }
    }
}

template<typename AOBasisSetType>
void AOQuartetIterator<AOBasisSetType>::increment_shells_() {
    auto old_shell_offsets = m_shell_quartet_iterator_.shell_offsets();
    auto old_shells        = m_shell_quartet_iterator_.current_quartet();

    ++m_shell_quartet_iterator_;

    if(!m_shell_quartet_iterator_) {
        return; // No more shell quartets, so we're done
    }

    auto new_shell_offsets = m_shell_quartet_iterator_.shell_offsets();

    for(size_type i = 0; i < 4; ++i) {
        auto new_i = new_shell_offsets[i];
        auto old_i = old_shell_offsets[i];

        if(new_i == old_i + 1) { // Shell i was incremented
            m_ao_offsets_[i] += old_shells[i].size();
        } else if(new_i == 0) {   // Shell i was reset
            m_ao_offsets_[i] = 0; // We wrapped around so reset
        } else if(new_i != old_i) {
            throw std::logic_error(
              "Expected shell to stay the same, increment by 1, or "
              "reset to 0, but got something else.");
        } // Else clause hit when shell offsets are same, so we do
          // nothing
    }
}

} // namespace integrals::utils
