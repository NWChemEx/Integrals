/*
 * Copyright 2024 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once
#include <simde/simde.hpp>

namespace integrals::ao_integrals {

class LibIntVisitor : public chemist::qm_operator::OperatorVisitor {
public:
    using s_e_type  = simde::type::s_e_type;
    using t_e_type  = simde::type::t_e_type;
    using v_ee_type = simde::type::v_ee_type;
    using v_en_type = simde::type::v_en_type;

    LibIntVisitor(const std::vector<libint2::BasisSet>& bases, double thresh,
                  std::size_t deriv = 0) :
      m_bases(bases), m_thresh(thresh), m_deriv(deriv) {};

    void run(const s_e_type& S_e) {
        m_engine = detail_::make_engine(m_bases, S_e, m_thresh, m_deriv);
    }

    void run(const t_e_type& T_e) {
        m_engine = detail_::make_engine(m_bases, T_e, m_thresh, m_deriv);
    }

    void run(const v_en_type& V_en) {
        m_engine = detail_::make_engine(m_bases, V_en, m_thresh, m_deriv);
    }

    void run(const v_ee_type& V_ee) {
        m_engine = detail_::make_engine(m_bases, V_ee, m_thresh, m_deriv);
    }

    libint2::Engine& engine() { return m_engine; }

private:
    const std::vector<libint2::BasisSet>& m_bases;
    double m_thresh;
    std::size_t m_deriv;
    libint2::Engine m_engine;
};

} // namespace integrals::ao_integrals