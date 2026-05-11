// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../registry.hpp"

#ifdef ATOMISTICA_HAVE_NETCDF
#include <netcdf.h>

namespace atomistica {

class OutputNC : public Hook {
    int         freq_;
    std::string file_;
    bool        write_velocities_;
    int         ncid_           = -1;
    int         frame_dim_      = -1;
    int         atom_dim_       = -1;
    int         spatial_dim_    = -1;
    int         time_var_       = -1;
    int         coords_var_     = -1;
    int         cell_lengths_var_ = -1;
    int         cell_angles_var_  = -1;
    int         vel_var_        = -1;
    size_t      frame_idx_      = 0;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 53; }
    std::string name() const override { return "OutputNC"; }

    explicit OutputNC(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 10))
        , file_(cfg.get_or<std::string>("file", "trajectory.nc"))
        , write_velocities_(cfg.get_or<bool>("velocities", false))
    {}

    ~OutputNC() {
        if (ncid_ >= 0)
            nc_close(ncid_);
    }

    void bind_to(AtomicSystem& sys, NeighborList&) override {
        check_nc(nc_create(file_.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid_),
                 "nc_create");

        // Dimensions
        check_nc(nc_def_dim(ncid_, "frame",   NC_UNLIMITED,                   &frame_dim_),   "def_dim frame");
        check_nc(nc_def_dim(ncid_, "atom",    static_cast<size_t>(sys.num_atoms()), &atom_dim_),    "def_dim atom");
        check_nc(nc_def_dim(ncid_, "spatial", 3,                              &spatial_dim_), "def_dim spatial");

        // time(frame)
        {
            int dims[1] = { frame_dim_ };
            check_nc(nc_def_var(ncid_, "time", NC_DOUBLE, 1, dims, &time_var_), "def_var time");
            const char* units = "picosecond";
            nc_put_att_text(ncid_, time_var_, "units", std::strlen(units), units);
        }

        // cell_lengths(frame, spatial)
        {
            int dims[2] = { frame_dim_, spatial_dim_ };
            check_nc(nc_def_var(ncid_, "cell_lengths", NC_DOUBLE, 2, dims, &cell_lengths_var_),
                     "def_var cell_lengths");
            const char* units = "angstrom";
            nc_put_att_text(ncid_, cell_lengths_var_, "units", std::strlen(units), units);
        }

        // cell_angles(frame, spatial)
        {
            int dims[2] = { frame_dim_, spatial_dim_ };
            check_nc(nc_def_var(ncid_, "cell_angles", NC_DOUBLE, 2, dims, &cell_angles_var_),
                     "def_var cell_angles");
            const char* units = "degree";
            nc_put_att_text(ncid_, cell_angles_var_, "units", std::strlen(units), units);
        }

        // coordinates(frame, atom, spatial)
        {
            int dims[3] = { frame_dim_, atom_dim_, spatial_dim_ };
            check_nc(nc_def_var(ncid_, "coordinates", NC_DOUBLE, 3, dims, &coords_var_),
                     "def_var coordinates");
            const char* units = "angstrom";
            nc_put_att_text(ncid_, coords_var_, "units", std::strlen(units), units);
        }

        // velocities(frame, atom, spatial)  — optional
        if (write_velocities_) {
            int dims[3] = { frame_dim_, atom_dim_, spatial_dim_ };
            check_nc(nc_def_var(ncid_, "velocities", NC_DOUBLE, 3, dims, &vel_var_),
                     "def_var velocities");
            const char* units = "angstrom/picosecond";
            nc_put_att_text(ncid_, vel_var_, "units", std::strlen(units), units);
        }

        // Global attributes (AMBER convention)
        const char* conv = "AMBER";
        nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", std::strlen(conv), conv);
        const char* ver = "1.0";
        nc_put_att_text(ncid_, NC_GLOBAL, "ConventionVersion", std::strlen(ver), ver);
        const char* prog = "Atomistica";
        nc_put_att_text(ncid_, NC_GLOBAL, "program", std::strlen(prog), prog);

        check_nc(nc_enddef(ncid_), "enddef");
    }

    void invoke(SimulationContext& ctx) override {
        if (ctx.step % freq_ != 0) return;

        size_t n = ctx.system.num_atoms();
        const Mat3& H = ctx.system.cell();

        // time
        {
            size_t start[1] = { frame_idx_ };
            size_t count[1] = { 1 };
            double t = ctx.time;
            check_nc(nc_put_vara_double(ncid_, time_var_, start, count, &t), "put time");
        }

        // cell_lengths: magnitudes of the three lattice vectors (columns of H)
        {
            size_t start[2] = { frame_idx_, 0 };
            size_t count[2] = { 1, 3 };
            double lengths[3] = {
                H.col(0).norm(),
                H.col(1).norm(),
                H.col(2).norm()
            };
            check_nc(nc_put_vara_double(ncid_, cell_lengths_var_, start, count, lengths),
                     "put cell_lengths");
        }

        // cell_angles: angles between lattice vectors (degrees)
        {
            size_t start[2] = { frame_idx_, 0 };
            size_t count[2] = { 1, 3 };
            Vec3 a = H.col(0), b = H.col(1), c = H.col(2);
            double la = a.norm(), lb = b.norm(), lc = c.norm();
            double alpha = (lb > 0.0 && lc > 0.0)
                ? std::acos(b.dot(c) / (lb * lc)) * (180.0 / M_PI) : 90.0;
            double beta  = (la > 0.0 && lc > 0.0)
                ? std::acos(a.dot(c) / (la * lc)) * (180.0 / M_PI) : 90.0;
            double gamma = (la > 0.0 && lb > 0.0)
                ? std::acos(a.dot(b) / (la * lb)) * (180.0 / M_PI) : 90.0;
            double angles[3] = { alpha, beta, gamma };
            check_nc(nc_put_vara_double(ncid_, cell_angles_var_, start, count, angles),
                     "put cell_angles");
        }

        // coordinates
        {
            std::vector<double> coords(n * 3);
            for (size_t i = 0; i < n; ++i) {
                Vec3 r = ctx.system.position(i);
                coords[3*i + 0] = r[0];
                coords[3*i + 1] = r[1];
                coords[3*i + 2] = r[2];
            }
            size_t start[3] = { frame_idx_, 0, 0 };
            size_t count[3] = { 1, n, 3 };
            check_nc(nc_put_vara_double(ncid_, coords_var_, start, count, coords.data()),
                     "put coordinates");
        }

        // velocities
        if (write_velocities_) {
            std::vector<double> vels(n * 3);
            for (size_t i = 0; i < n; ++i) {
                Vec3 v = ctx.system.velocity(i);
                vels[3*i + 0] = v[0];
                vels[3*i + 1] = v[1];
                vels[3*i + 2] = v[2];
            }
            size_t start[3] = { frame_idx_, 0, 0 };
            size_t count[3] = { 1, n, 3 };
            check_nc(nc_put_vara_double(ncid_, vel_var_, start, count, vels.data()),
                     "put velocities");
        }

        nc_sync(ncid_);
        ++frame_idx_;
    }

private:
    void check_nc(int status, const char* msg) {
        if (status != NC_NOERR)
            throw std::runtime_error(
                std::string("OutputNC: ") + msg + ": " + nc_strerror(status));
    }
};

REGISTER_HOOK("OutputNC", OutputNC)

} // namespace atomistica

#else // ATOMISTICA_HAVE_NETCDF

namespace atomistica {

class OutputNC : public Hook {
public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 53; }
    std::string name() const override { return "OutputNC"; }

    explicit OutputNC(const Config&) {
        throw std::runtime_error(
            "OutputNC requires NetCDF (rebuild with NetCDF enabled)");
    }

    void invoke(SimulationContext&) override {}
};

REGISTER_HOOK("OutputNC", OutputNC)

} // namespace atomistica

#endif // ATOMISTICA_HAVE_NETCDF
