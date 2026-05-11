// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
//
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ======================================================================

#pragma once

#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../config.hpp"
#include "../../../include/atomistica/config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../registry.hpp"
#include "../../atoms_io.hpp"   // for Z_to_symbol / symbol_to_Z

namespace atomistica {

// Charge equilibration (QEq / EEM) via the electronegativity equalization
// method (Rappe & Goddard 1991).
//
// Energy model:
//   E = sum_i [chi_i * q_i + 0.5 * J_i * q_i^2] + E_Coulomb(q)
//
// Equilibrium condition (with Lagrange multiplier mu for charge neutrality):
//   chi_i + J_i * q_i + phi_i = mu   for all i
//
// Self-consistent solution (Anderson linear mixing):
//   q_i^new = (mu - chi_i - phi_i) / J_i
//   q       = (1-mixing)*q + mixing*q_new
//
// Default QEq parameters (Rappe & Goddard 1991, Table 1):
//   H:  chi=4.528 eV/e,  J=13.890 eV/e^2
//   C:  chi=5.343 eV/e,  J=10.126 eV/e^2
//   N:  chi=6.899 eV/e,  J=11.760 eV/e^2
//   O:  chi=8.741 eV/e,  J=13.364 eV/e^2
//   Si: chi=4.168 eV/e,  J= 6.974 eV/e^2
//   S:  chi=6.921 eV/e,  J= 8.972 eV/e^2
//   P:  chi=5.463 eV/e,  J= 8.000 eV/e^2
//   Cl: chi=8.564 eV/e,  J= 9.892 eV/e^2
//   Fe: chi=4.060 eV/e,  J= 8.045 eV/e^2
//   Cu: chi=4.200 eV/e,  J= 6.997 eV/e^2
//
// Config keys:
//   convergence    — charge change tolerance (default 1e-4 e)
//   max_iter       — maximum SCF iterations (default 200)
//   mixing         — linear mixing parameter (default 0.2)
//   total_charge   — constraint: sum(q_i) = total_charge (default 0)
//   chi_<ELEM>     — per-element electronegativity override (eV/e)
//   J_<ELEM>       — per-element hardness override (eV/e^2)
//
// The Coulomb solver must be registered in the same simulation (in the
// SimulationContext::coulomb slot).  If no Coulomb solver is present,
// this hook is a no-op.
class ChargeEquilibration : public Hook {
    double convergence_  = 1e-4;
    int    max_iter_     = 200;
    double mixing_       = 0.2;
    double total_charge_ = 0.0;

    // Per-element QEq parameters (keyed by atomic number Z)
    std::map<int, double> chi_;
    std::map<int, double> J_;

    // Working arrays (allocated in bind_to)
    ArrayX q_;
    ArrayX phi_;

    // Default Rappe-Goddard parameters
    static std::map<int, double> default_chi() {
        return {{1, 4.528}, {6, 5.343}, {7, 6.899}, {8, 8.741},
                {14, 4.168}, {15, 5.463}, {16, 6.921}, {17, 8.564},
                {26, 4.060}, {29, 4.200}};
    }
    static std::map<int, double> default_J() {
        return {{1, 13.890}, {6, 10.126}, {7, 11.760}, {8, 13.364},
                {14, 6.974}, {15, 8.000}, {16, 8.972}, {17, 9.892},
                {26, 8.045}, {29, 6.997}};
    }

public:
    HookPoint   hook_point() const override { return HookPoint::PRE_FORCE; }
    int         priority()   const override { return 5; }
    std::string name()       const override { return "ChargeEquilibration"; }

    explicit ChargeEquilibration(const Config& cfg)
        : convergence_ (cfg.get_or<double>("convergence",   1e-4))
        , max_iter_    (cfg.get_or<int>   ("max_iter",       200))
        , mixing_      (cfg.get_or<double>("mixing",         0.2))
        , total_charge_(cfg.get_or<double>("total_charge",   0.0))
        , chi_(default_chi())
        , J_  (default_J())
    {
        // Read per-element overrides: chi_C, J_C, chi_H, J_H, …
        static const char* symbols[] = {
            "H","He","Li","Be","B","C","N","O","F","Ne",
            "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
            "Fe","Ni","Cu","Ag","W","Pt","Au", nullptr
        };
        for (int s = 0; symbols[s]; ++s) {
            std::string sym(symbols[s]);
            std::string chi_key = "chi_" + sym;
            std::string J_key   = "J_"   + sym;
            if (cfg.has_key(chi_key)) {
                int Z = symbol_to_Z(sym);
                chi_[Z] = cfg.get_or<double>(chi_key, chi_[Z]);
            }
            if (cfg.has_key(J_key)) {
                int Z = symbol_to_Z(sym);
                J_[Z] = cfg.get_or<double>(J_key, J_[Z]);
            }
        }
    }

    void bind_to(AtomicSystem& sys, NeighborList&) override {
        size_t n = sys.num_atoms();
        q_.resize(n);
        phi_.resize(n);
        q_.setZero();

        if (!sys.properties().has("charges"))
            sys.properties().add<ArrayX>("charges", n);
    }

    void invoke(SimulationContext& ctx) override {
        if (!ctx.coulomb) return;

        int n = static_cast<int>(ctx.system.num_atoms());

        // Warm-start from previously converged charges if available.
        if (ctx.system.properties().has("charges")) {
            q_ = ctx.system.properties().get<ArrayX>("charges");
            if (static_cast<int>(q_.size()) != n) q_.resize(n), q_.setZero();
        }

        // SCF loop
        for (int iter = 0; iter < max_iter_; ++iter) {
            ctx.coulomb->compute_potential(ctx.system, ctx.nl, q_, phi_);

            // Lagrange multiplier enforcing sum(q_i) = total_charge_
            double sum_inv_J = 0.0;
            double sum_num   = 0.0;
            for (int i = 0; i < n; ++i) {
                int    Z  = ctx.system.atomic_number(static_cast<size_t>(i));
                double Ji = J_.count(Z) ? J_.at(Z) : 10.0;
                double xi = chi_.count(Z) ? chi_.at(Z) : 5.0;
                sum_inv_J += 1.0 / Ji;
                sum_num   += (xi + static_cast<double>(phi_[i])) / Ji;
            }
            if (sum_inv_J < 1e-30) break;
            double mu = (total_charge_ + sum_num) / sum_inv_J;

            // New charges and convergence check
            double dq_max = 0.0;
            ArrayX q_new(n);
            for (int i = 0; i < n; ++i) {
                int    Z  = ctx.system.atomic_number(static_cast<size_t>(i));
                double Ji = J_.count(Z) ? J_.at(Z) : 10.0;
                double xi = chi_.count(Z) ? chi_.at(Z) : 5.0;
                q_new[i] = static_cast<Scalar>(
                    (mu - xi - static_cast<double>(phi_[i])) / Ji);
                dq_max = std::max(dq_max, std::abs(
                    static_cast<double>(q_new[i] - q_[i])));
            }

            q_ = (1.0 - mixing_) * q_ + mixing_ * q_new;

            if (dq_max < convergence_) break;
        }

        ctx.system.properties().get<ArrayX>("charges") = q_;
    }
};

REGISTER_HOOK("ChargeEquilibration", ChargeEquilibration)

} // namespace atomistica
