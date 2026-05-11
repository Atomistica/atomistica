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

#include <string>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../simulation_context.hpp"
#include "../hook.hpp"
#include "../integrator.hpp"
#include "../md_utils.hpp"

namespace atomistica {

class VelocityVerlet : public Integrator {
public:
    explicit VelocityVerlet(const Config&) {}

    std::string name() const override { return "VelocityVerlet"; }

    bool step(SimulationContext& ctx) override {
        const double dt = ctx.dt;
        const size_t n = ctx.system.num_atoms();

        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3 f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            Vec3 v = ctx.system.velocity(i);
            ctx.system.set_velocity(i, v + 0.5 * dt * f / m);
        }

        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            Vec3 r = ctx.system.position(i);
            Vec3 v = ctx.system.velocity(i);
            ctx.system.set_position(i, r + dt * v);
        }

        wrap_positions(ctx.system);

        ctx.invoke_hooks(HookPoint::PRE_FORCE);
        ctx.nl.update(ctx.system);
        ctx.compute_forces();
        ctx.invoke_hooks(HookPoint::POST_FORCE);

        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3 f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            Vec3 v = ctx.system.velocity(i);
            ctx.system.set_velocity(i, v + 0.5 * dt * f / m);
        }

        ctx.invoke_hooks(HookPoint::POST_STEP);

        ctx.time += ctx.dt;
        ++ctx.step;

        return true;
    }
};

} // namespace atomistica
