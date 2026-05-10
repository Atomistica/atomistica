// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"

#include "../../registry.hpp"

namespace atomistica {

// ---------------------------------------------------------------------------
// FreezeAtoms: zero forces and velocities of frozen atoms (group <= 0).
// ---------------------------------------------------------------------------

class FreezeAtoms : public Hook {
public:
    HookPoint hook_point() const override { return HookPoint::POST_FORCE; }
    int priority() const override { return 9; }
    std::string name() const override { return "FreezeAtoms"; }

    explicit FreezeAtoms(const Config&) {}

    void invoke(SimulationContext& ctx) override {
        for (size_t i = 0; i < ctx.system.num_atoms(); ++i) {
            if (is_frozen(ctx.system, i)) {
                ctx.system.forces().col(i).setZero();
                ctx.system.set_velocity(i, Vec3::Zero());
            }
        }
    }
};

REGISTER_HOOK("FreezeAtoms", FreezeAtoms)

// ---------------------------------------------------------------------------
// SETTLE: SHAKE-based rigid-water constraint.
// ---------------------------------------------------------------------------

struct WaterMol { int O, H1, H2; };

class SETTLE : public Hook {
    double d_OH_;
    double d_HH_;
    double tol_;
    int    max_iter_;
    std::vector<WaterMol> molecules_;

public:
    HookPoint hook_point() const override { return HookPoint::PRE_FORCE; }
    int priority() const override { return 15; }
    std::string name() const override { return "SETTLE"; }

    explicit SETTLE(const Config& cfg)
        : d_OH_(cfg.get_or<double>("d_OH", 0.9572))
        , d_HH_(cfg.get_or<double>("d_HH", 1.5139))
        , tol_(cfg.get_or<double>("tol", 1e-6))
        , max_iter_(cfg.get_or<int>("max_iter", 100))
    {}

    void bind_to(AtomicSystem& sys, NeighborList&) override {
        molecules_.clear();
        size_t n = sys.num_atoms();
        double search_r2 = 1.2 * 1.2;

        for (size_t i = 0; i < n; ++i) {
            if (sys.atomic_number(i) != 8) continue;
            Vec3 rO = sys.position(i);

            // Find up to two nearest H atoms within 1.2 Å
            int h1 = -1, h2 = -1;
            double d1 = search_r2, d2 = search_r2;

            for (size_t j = 0; j < n; ++j) {
                if (sys.atomic_number(j) != 1) continue;
                Vec3 dr = sys.position(j) - rO;
                // Use minimum image for PBC
                dr = sys.minimum_image(dr);
                double d2j = dr.squaredNorm();
                if (d2j < d1) {
                    d2  = d1; h2  = h1;
                    d1  = d2j; h1 = static_cast<int>(j);
                } else if (d2j < d2) {
                    d2 = d2j; h2 = static_cast<int>(j);
                }
            }

            if (h1 >= 0 && h2 >= 0)
                molecules_.push_back({static_cast<int>(i), h1, h2});
        }

        std::cout << "SETTLE: found " << molecules_.size()
                  << " water molecule(s).\n";
    }

    void invoke(SimulationContext& ctx) override {
        double d_OH2 = d_OH_ * d_OH_;
        double d_HH2 = d_HH_ * d_HH_;

        for (const auto& mol : molecules_) {
            double mO  = ctx.system.mass(static_cast<size_t>(mol.O));
            double mH1 = ctx.system.mass(static_cast<size_t>(mol.H1));
            double mH2 = ctx.system.mass(static_cast<size_t>(mol.H2));

            double inv_mO  = 1.0 / mO;
            double inv_mH1 = 1.0 / mH1;
            double inv_mH2 = 1.0 / mH2;

            for (int iter = 0; iter < max_iter_; ++iter) {
                Vec3 rO  = ctx.system.position(static_cast<size_t>(mol.O));
                Vec3 rH1 = ctx.system.position(static_cast<size_t>(mol.H1));
                Vec3 rH2 = ctx.system.position(static_cast<size_t>(mol.H2));

                // Constraint 1: O-H1
                Vec3   dr_OH1 = rO - rH1;
                double d2_OH1 = dr_OH1.squaredNorm();
                double lam1   = (d_OH2 - d2_OH1)
                                / (2.0 * (inv_mO + inv_mH1) * d2_OH1);
                Vec3 drO_1  =  lam1 * inv_mO  * dr_OH1;
                Vec3 drH1_1 = -lam1 * inv_mH1 * dr_OH1;

                // Constraint 2: O-H2
                Vec3   dr_OH2 = rO - rH2;
                double d2_OH2 = dr_OH2.squaredNorm();
                double lam2   = (d_OH2 - d2_OH2)
                                / (2.0 * (inv_mO + inv_mH2) * d2_OH2);
                Vec3 drO_2  =  lam2 * inv_mO  * dr_OH2;
                Vec3 drH2_2 = -lam2 * inv_mH2 * dr_OH2;

                // Constraint 3: H1-H2
                Vec3   dr_HH = rH1 - rH2;
                double d2_HH = dr_HH.squaredNorm();
                double lam3  = (d_HH2 - d2_HH)
                               / (2.0 * (inv_mH1 + inv_mH2) * d2_HH);
                Vec3 drH1_3 =  lam3 * inv_mH1 * dr_HH;
                Vec3 drH2_3 = -lam3 * inv_mH2 * dr_HH;

                ctx.system.set_position(static_cast<size_t>(mol.O),
                    rO  + drO_1  + drO_2);
                ctx.system.set_position(static_cast<size_t>(mol.H1),
                    rH1 + drH1_1 + drH1_3);
                ctx.system.set_position(static_cast<size_t>(mol.H2),
                    rH2 + drH2_2 + drH2_3);

                double max_delta = std::max({std::abs(lam1), std::abs(lam2), std::abs(lam3)});
                if (max_delta < tol_) break;
            }
        }

        if (!molecules_.empty())
            ctx.system.positions_changed();
    }
};

REGISTER_HOOK("SETTLE", SETTLE)

} // namespace atomistica
