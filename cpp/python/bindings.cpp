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

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <atomistica/atomistica.hpp>
#include <atomistica/potentials/eam/eam.hpp>

namespace py = pybind11;
using namespace atomistica;

PYBIND11_MODULE(_atomistica_cpp, m) {
    m.doc() = "Atomistica C++ - Interatomic potentials library";

    // PotentialResults
    py::class_<PotentialResults>(m, "PotentialResults")
        .def(py::init<>())
        .def_readwrite("energy", &PotentialResults::energy)
        .def_readwrite("virial", &PotentialResults::virial)
        .def_property_readonly("has_energy_per_atom",
            [](const PotentialResults& r) { return r.energy_per_atom.has_value(); })
        .def_property_readonly("energy_per_atom",
            [](const PotentialResults& r) -> py::object {
                if (r.energy_per_atom.has_value()) {
                    return py::cast(*r.energy_per_atom);
                }
                return py::none();
            });

    // AtomicSystem
    py::class_<AtomicSystem>(m, "AtomicSystem")
        .def(py::init<>())
        .def(py::init<std::size_t>(), py::arg("num_atoms"))
        .def("resize", &AtomicSystem::resize)
        .def_property_readonly("num_atoms", &AtomicSystem::num_atoms)
        .def_property("cell",
            [](const AtomicSystem& s) { return s.cell(); },
            [](AtomicSystem& s, const Mat3& c) { s.set_cell(c); })
        .def_property("pbc",
            [](const AtomicSystem& s) { return s.pbc(); },
            [](AtomicSystem& s, const std::array<bool, 3>& p) { s.pbc() = p; })
        .def_property("positions",
            [](AtomicSystem& s) -> Eigen::Ref<Array3X> { return s.positions(); },
            [](AtomicSystem& s, const Array3X& p) { s.positions() = p; s.positions_changed(); })
        .def_property("atomic_numbers",
            [](AtomicSystem& s) -> Eigen::Ref<ArrayXi> { return s.atomic_numbers(); },
            [](AtomicSystem& s, const ArrayXi& z) { s.atomic_numbers() = z; })
        .def_property_readonly("forces",
            [](AtomicSystem& s) -> Eigen::Ref<Array3X> { return s.forces(); })
        .def("zero_forces", &AtomicSystem::zero_forces)
        .def("volume", &AtomicSystem::volume)
        .def("minimum_image", &AtomicSystem::minimum_image)
        .def("wrap_position", &AtomicSystem::wrap_position)
        .def("positions_changed", &AtomicSystem::positions_changed)
        .def("cell_changed", &AtomicSystem::cell_changed);

    // Neighbor
    py::class_<Neighbor>(m, "Neighbor")
        .def_readonly("index", &Neighbor::index)
        .def_readonly("cell_shift", &Neighbor::cell_shift);

    // NeighborList
    py::class_<NeighborList>(m, "NeighborList")
        .def(py::init<>())
        .def("set_cutoff", &NeighborList::set_cutoff)
        .def_property_readonly("cutoff", &NeighborList::cutoff)
        .def("set_verlet_shell", &NeighborList::set_verlet_shell)
        .def_property_readonly("verlet_shell", &NeighborList::verlet_shell)
        .def("update", &NeighborList::update)
        .def("invalidate", &NeighborList::invalidate)
        .def_property_readonly("num_atoms", &NeighborList::num_atoms)
        .def("num_neighbors", &NeighborList::num_neighbors)
        .def_property_readonly("num_pairs", &NeighborList::num_pairs)
        .def("neighbors", [](const NeighborList& nl, std::size_t i) -> std::vector<Neighbor> {
            auto [begin, end] = nl.neighbors(i);
            return std::vector<Neighbor>(begin, end);
        });

    // LJCut potential
    py::class_<LJCut>(m, "LJCut")
        .def(py::init<>())
        .def(py::init<int, Scalar, Scalar, Scalar>(),
             py::arg("Z"), py::arg("epsilon"), py::arg("sigma"), py::arg("cutoff"))
        .def("set_params", &LJCut::set_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("epsilon"),
             py::arg("sigma"), py::arg("cutoff"))
        .def("cutoff", &LJCut::cutoff)
        .def("bind_to", &LJCut::bind_to)
        .def("compute", &LJCut::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true);

    // LJCutShift potential
    py::class_<LJCutShift>(m, "LJCutShift")
        .def(py::init<>())
        .def(py::init<int, Scalar, Scalar, Scalar>(),
             py::arg("Z"), py::arg("epsilon"), py::arg("sigma"), py::arg("cutoff"))
        .def("set_params", &LJCutShift::set_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("epsilon"),
             py::arg("sigma"), py::arg("cutoff"))
        .def("cutoff", &LJCutShift::cutoff)
        .def("bind_to", &LJCutShift::bind_to)
        .def("compute", &LJCutShift::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true);

    // Spline classes
    py::class_<CubicSpline>(m, "CubicSpline")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar, const std::vector<Scalar>&>(),
             py::arg("x_min"), py::arg("x_max"), py::arg("y"))
        .def("init", &CubicSpline::init)
        .def("is_valid", &CubicSpline::is_valid)
        .def("eval", [](const CubicSpline& s, Scalar x) {
            auto r = s.eval(x);
            return std::make_pair(r.value, r.derivative);
        })
        .def("value", &CubicSpline::value)
        .def_property_readonly("x_min", &CubicSpline::x_min)
        .def_property_readonly("x_max", &CubicSpline::x_max);

    py::class_<NonUniformSpline>(m, "NonUniformSpline")
        .def(py::init<>())
        .def(py::init<const std::vector<Scalar>&, const std::vector<Scalar>&>(),
             py::arg("x"), py::arg("y"))
        .def("init", &NonUniformSpline::init)
        .def("is_valid", &NonUniformSpline::is_valid)
        .def("eval", [](const NonUniformSpline& s, Scalar x) {
            auto r = s.eval(x);
            return std::make_pair(r.value, r.derivative);
        })
        .def("value", &NonUniformSpline::value)
        .def_property_readonly("x_min", &NonUniformSpline::x_min)
        .def_property_readonly("x_max", &NonUniformSpline::x_max);

    // =========================================================================
    // Cutoff functions
    // =========================================================================

    // CutoffResult helper
    py::class_<CutoffResult>(m, "CutoffResult")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar>(), py::arg("fc"), py::arg("dfc"))
        .def_readwrite("fc", &CutoffResult::fc)
        .def_readwrite("dfc", &CutoffResult::dfc);

    // TrigOffCutoff
    py::class_<TrigOffCutoff>(m, "TrigOffCutoff")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar>(), py::arg("r1"), py::arg("r2"))
        .def("init", &TrigOffCutoff::init)
        .def_property_readonly("r1", &TrigOffCutoff::r1)
        .def_property_readonly("r2", &TrigOffCutoff::r2)
        .def_property_readonly("cutoff", &TrigOffCutoff::cutoff)
        .def("__call__", &TrigOffCutoff::operator());

    // TrigOnCutoff
    py::class_<TrigOnCutoff>(m, "TrigOnCutoff")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar>(), py::arg("r1"), py::arg("r2"))
        .def("init", &TrigOnCutoff::init)
        .def_property_readonly("r1", &TrigOnCutoff::r1)
        .def_property_readonly("r2", &TrigOnCutoff::r2)
        .def("__call__", &TrigOnCutoff::operator());

    // ExpCutoff
    py::class_<ExpCutoff>(m, "ExpCutoff")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar>(), py::arg("r1"), py::arg("r2"))
        .def("init", &ExpCutoff::init)
        .def_property_readonly("r1", &ExpCutoff::r1)
        .def_property_readonly("r2", &ExpCutoff::r2)
        .def_property_readonly("cutoff", &ExpCutoff::cutoff)
        .def("__call__", &ExpCutoff::operator());

    // =========================================================================
    // BOP Potentials
    // =========================================================================

    // TersoffElementParams
    py::class_<TersoffElementParams>(m, "TersoffElementParams")
        .def(py::init<>())
        .def_readwrite("beta", &TersoffElementParams::beta)
        .def_readwrite("n", &TersoffElementParams::n)
        .def_readwrite("c", &TersoffElementParams::c)
        .def_readwrite("d", &TersoffElementParams::d)
        .def_readwrite("h", &TersoffElementParams::h)
        .def("precompute_angular", &TersoffElementParams::precompute_angular);

    // TersoffPairParams
    py::class_<TersoffPairParams>(m, "TersoffPairParams")
        .def(py::init<>())
        .def_readwrite("A", &TersoffPairParams::A)
        .def_readwrite("B", &TersoffPairParams::B)
        .def_readwrite("lambda_", &TersoffPairParams::lambda)
        .def_readwrite("mu", &TersoffPairParams::mu)
        .def_readwrite("r1", &TersoffPairParams::r1)
        .def_readwrite("r2", &TersoffPairParams::r2)
        .def_readwrite("chi", &TersoffPairParams::chi)
        .def("init_cutoff", &TersoffPairParams::init_cutoff);

    // Tersoff potential (non-screened)
    py::class_<Tersoff<false>>(m, "Tersoff")
        .def(py::init<>())
        .def("add_element", &Tersoff<false>::add_element,
             py::arg("Z"), py::arg("params"),
             "Add element with given atomic number and parameters")
        .def("set_pair_params", &Tersoff<false>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"),
             "Set pair parameters for element pair")
        .def("load_parameters", &Tersoff<false>::load_parameters,
             py::arg("name"),
             "Load built-in parameter set by name")
        .def("cutoff", &Tersoff<false>::cutoff,
             "Get maximum cutoff radius")
        .def("num_elements", &Tersoff<false>::num_elements,
             "Get number of elements defined")
        .def("element_index", &Tersoff<false>::element_index,
             py::arg("Z"),
             "Get internal element index for atomic number Z (-1 if not found)")
        .def("pair_type", &Tersoff<false>::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type index for element pair")
        .def("compute", &Tersoff<false>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        // Expose internal functions for testing/analysis
        .def("repulsive", [](const Tersoff<false>& pot, int ptype, Scalar r) {
            auto [V, dV] = pot.repulsive(ptype, r);
            return std::make_pair(V, dV);
        }, py::arg("ptype"), py::arg("r"))
        .def("attractive", [](const Tersoff<false>& pot, int ptype, Scalar r) {
            auto [V, dV] = pot.attractive(ptype, r);
            return std::make_pair(V, dV);
        }, py::arg("ptype"), py::arg("r"))
        .def("angular_function", [](const Tersoff<false>& pot,
                int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) {
            auto [g, dg] = pot.angular_function(eli, elj, elk, ptype_ij, ptype_ik, cos_theta);
            return std::make_pair(g, dg);
        }, py::arg("eli"), py::arg("elj"), py::arg("elk"),
           py::arg("ptype_ij"), py::arg("ptype_ik"), py::arg("cos_theta"))
        .def("bond_order", [](const Tersoff<false>& pot, int eli, int ptype, Scalar z) {
            auto [b, db] = pot.bond_order(eli, ptype, z);
            return std::make_pair(b, db);
        }, py::arg("eli"), py::arg("ptype"), py::arg("z"));

    // ScreeningParams
    py::class_<ScreeningParams>(m, "ScreeningParams")
        .def(py::init<>())
        .def_readwrite("Cmin", &ScreeningParams::Cmin)
        .def_readwrite("Cmax", &ScreeningParams::Cmax)
        .def_readwrite("cut_in_l", &ScreeningParams::cut_in_l)
        .def_readwrite("cut_in_h", &ScreeningParams::cut_in_h)
        .def_readwrite("cut_out_l", &ScreeningParams::cut_out_l)
        .def_readwrite("cut_out_h", &ScreeningParams::cut_out_h)
        .def_readwrite("cut_bo_l", &ScreeningParams::cut_bo_l)
        .def_readwrite("cut_bo_h", &ScreeningParams::cut_bo_h)
        .def("precompute", &ScreeningParams::precompute);

    // Screened Tersoff potential
    py::class_<Tersoff<true>>(m, "TersoffScr")
        .def(py::init<>())
        .def("add_element", &Tersoff<true>::add_element,
             py::arg("Z"), py::arg("params"),
             "Add element with given atomic number and parameters")
        .def("set_pair_params", &Tersoff<true>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"),
             "Set pair parameters for element pair")
        .def("load_parameters", &Tersoff<true>::load_parameters,
             py::arg("name"),
             "Load built-in parameter set by name")
        .def("cutoff", &Tersoff<true>::cutoff,
             "Get maximum cutoff radius")
        .def("num_elements", &Tersoff<true>::num_elements,
             "Get number of elements defined")
        .def("element_index", &Tersoff<true>::element_index,
             py::arg("Z"),
             "Get internal element index for atomic number Z (-1 if not found)")
        .def("pair_type", &Tersoff<true>::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type index for element pair")
        .def("compute", &Tersoff<true>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        .def("screening_params", &Tersoff<true>::screening_params,
             py::arg("ptype"),
             "Get screening parameters for pair type",
             py::return_value_policy::reference_internal);

    // Available parameter sets
    m.def("available_tersoff_parameters", []() {
        return std::vector<std::string>{"Tersoff_PRB_39_5566_Si_C"};
    }, "List available built-in Tersoff parameter sets");

    // =========================================================================
    // EAM Potentials
    // =========================================================================

    // EAMElementInfo
    py::class_<EAMElementInfo>(m, "EAMElementInfo")
        .def(py::init<>())
        .def_readwrite("symbol", &EAMElementInfo::symbol)
        .def_readwrite("atomic_number", &EAMElementInfo::atomic_number)
        .def_readwrite("mass", &EAMElementInfo::mass)
        .def_readwrite("lattice_constant", &EAMElementInfo::lattice_constant)
        .def_readwrite("lattice_type", &EAMElementInfo::lattice_type);

    // TabulatedEAM (single element, funcfl format)
    py::class_<TabulatedEAM>(m, "TabulatedEAM")
        .def(py::init<>())
        .def("load", &TabulatedEAM::load,
             py::arg("filename"),
             "Load EAM potential from funcfl format file")
        .def("is_valid", &TabulatedEAM::is_valid,
             "Check if potential is loaded")
        .def("element_info", &TabulatedEAM::element_info,
             "Get element information",
             py::return_value_policy::reference_internal)
        .def("cutoff", &TabulatedEAM::cutoff,
             "Get cutoff radius")
        .def("embedding", [](const TabulatedEAM& pot, Scalar rho) {
            auto r = pot.embedding(rho);
            return std::make_pair(r.value, r.derivative);
        }, py::arg("rho"), "Evaluate embedding function F(rho)")
        .def("effective_charge", [](const TabulatedEAM& pot, Scalar r) {
            auto res = pot.effective_charge(r);
            return std::make_pair(res.value, res.derivative);
        }, py::arg("r"), "Evaluate effective charge Z(r)")
        .def("density", [](const TabulatedEAM& pot, Scalar r) {
            auto res = pot.density(r);
            return std::make_pair(res.value, res.derivative);
        }, py::arg("r"), "Evaluate electron density rho(r)")
        .def("pair_potential", [](const TabulatedEAM& pot, Scalar r) {
            auto res = pot.pair_potential(r);
            return std::make_pair(res.value, res.derivative);
        }, py::arg("r"), "Evaluate pair potential phi(r) = Z(r)^2 / r")
        .def("compute", &TabulatedEAM::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial");

    // TabulatedAlloyEAM (multi-element, setfl format)
    py::class_<TabulatedAlloyEAM>(m, "TabulatedAlloyEAM")
        .def(py::init<>())
        .def("load", &TabulatedAlloyEAM::load,
             py::arg("filename"),
             "Load EAM potential from setfl/alloy format file")
        .def("is_valid", &TabulatedAlloyEAM::is_valid,
             "Check if potential is loaded")
        .def("num_elements", &TabulatedAlloyEAM::num_elements,
             "Get number of elements")
        .def("element_info", &TabulatedAlloyEAM::element_info,
             py::arg("elem_idx"),
             "Get element information for element index",
             py::return_value_policy::reference_internal)
        .def("element_index", &TabulatedAlloyEAM::element_index,
             py::arg("symbol"),
             "Get element index from symbol (-1 if not found)")
        .def("element_index_by_Z", &TabulatedAlloyEAM::element_index_by_Z,
             py::arg("Z"),
             "Get element index from atomic number (-1 if not found)")
        .def("element_symbols", &TabulatedAlloyEAM::element_symbols,
             "Get list of element symbols")
        .def("cutoff", &TabulatedAlloyEAM::cutoff,
             "Get cutoff radius")
        .def("embedding", [](const TabulatedAlloyEAM& pot, int elem_idx, Scalar rho) {
            auto r = pot.embedding(elem_idx, rho);
            return std::make_pair(r.value, r.derivative);
        }, py::arg("elem_idx"), py::arg("rho"),
           "Evaluate embedding function F_i(rho) for element i")
        .def("density", [](const TabulatedAlloyEAM& pot, int elem_idx, Scalar r) {
            auto res = pot.density(elem_idx, r);
            return std::make_pair(res.value, res.derivative);
        }, py::arg("elem_idx"), py::arg("r"),
           "Evaluate electron density rho_i(r) for element i")
        .def("pair_potential", [](const TabulatedAlloyEAM& pot, int elem_i, int elem_j, Scalar r) {
            auto res = pot.pair_potential(elem_i, elem_j, r);
            return std::make_pair(res.value, res.derivative);
        }, py::arg("elem_i"), py::arg("elem_j"), py::arg("r"),
           "Evaluate pair potential phi_ij(r) between elements i and j")
        .def("compute", &TabulatedAlloyEAM::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial");

    // =========================================================================
    // Coulomb Potentials
    // =========================================================================

    // Coulomb constant
    m.attr("COULOMB_CONST") = COULOMB_CONST;

    // DirectCoulomb (O(N^2) for non-periodic systems)
    py::class_<DirectCoulomb>(m, "DirectCoulomb")
        .def(py::init<>())
        .def(py::init<Scalar>(), py::arg("epsilon_r") = 1.0,
             "Construct with optional dielectric constant")
        .def("set_epsilon_r", &DirectCoulomb::set_epsilon_r,
             py::arg("epsilon_r"),
             "Set relative dielectric constant")
        .def_property_readonly("epsilon_r", &DirectCoulomb::epsilon_r,
             "Get relative dielectric constant")
        .def("set_charges", &DirectCoulomb::set_charges,
             py::arg("charges"),
             "Set charges for all atoms (in units of e)")
        .def_property_readonly("charges", &DirectCoulomb::charges,
             "Get charges")
        .def("cutoff", &DirectCoulomb::cutoff,
             "Get cutoff (infinity for DirectCoulomb)")
        .def("compute", &DirectCoulomb::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial");

    // CutoffCoulomb (hard cutoff - use with care!)
    py::class_<CutoffCoulomb>(m, "CutoffCoulomb")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar>(),
             py::arg("cutoff") = 10.0, py::arg("epsilon_r") = 1.0,
             "Construct with cutoff and optional dielectric constant")
        .def("set_cutoff", &CutoffCoulomb::set_cutoff,
             py::arg("cutoff"),
             "Set cutoff radius")
        .def("set_epsilon_r", &CutoffCoulomb::set_epsilon_r,
             py::arg("epsilon_r"),
             "Set relative dielectric constant")
        .def_property_readonly("epsilon_r", &CutoffCoulomb::epsilon_r,
             "Get relative dielectric constant")
        .def("set_charges", &CutoffCoulomb::set_charges,
             py::arg("charges"),
             "Set charges for all atoms (in units of e)")
        .def_property_readonly("charges", &CutoffCoulomb::charges,
             "Get charges")
        .def("cutoff", &CutoffCoulomb::cutoff,
             "Get cutoff radius")
        .def("compute", &CutoffCoulomb::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial");

    // WolfCoulomb (damped shifted force method)
    py::class_<WolfCoulomb>(m, "WolfCoulomb")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar, Scalar>(),
             py::arg("cutoff") = 10.0, py::arg("alpha") = 0.0, py::arg("epsilon_r") = 1.0,
             "Construct with cutoff, damping parameter (0=auto), and dielectric constant")
        .def("set_cutoff", &WolfCoulomb::set_cutoff,
             py::arg("cutoff"),
             "Set cutoff radius")
        .def("set_alpha", &WolfCoulomb::set_alpha,
             py::arg("alpha"),
             "Set damping parameter (0 for auto-compute)")
        .def("set_epsilon_r", &WolfCoulomb::set_epsilon_r,
             py::arg("epsilon_r"),
             "Set relative dielectric constant")
        .def_property_readonly("alpha", &WolfCoulomb::alpha,
             "Get damping parameter")
        .def_property_readonly("epsilon_r", &WolfCoulomb::epsilon_r,
             "Get relative dielectric constant")
        .def("set_charges", &WolfCoulomb::set_charges,
             py::arg("charges"),
             "Set charges for all atoms (in units of e)")
        .def_property_readonly("charges", &WolfCoulomb::charges,
             "Get charges")
        .def("cutoff", &WolfCoulomb::cutoff,
             "Get cutoff radius")
        .def("compute", &WolfCoulomb::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial");

    // PMECoulomb (Particle Mesh Ewald for periodic systems)
    py::class_<PMECoulomb>(m, "PMECoulomb")
        .def(py::init<Scalar, int, int, int, int, Scalar>(),
             py::arg("cutoff") = 10.0,
             py::arg("grid_x") = 32, py::arg("grid_y") = 32, py::arg("grid_z") = 32,
             py::arg("order") = 4, py::arg("alpha") = 0.0,
             "Construct PME solver with cutoff, grid dimensions, B-spline order, and alpha (0=auto)")
        .def("set_cutoff", &PMECoulomb::set_cutoff,
             py::arg("cutoff"),
             "Set real-space cutoff")
        .def("set_alpha", &PMECoulomb::set_alpha,
             py::arg("alpha"),
             "Set Ewald parameter (0 for auto-compute)")
        .def("set_grid", &PMECoulomb::set_grid,
             py::arg("grid_x"), py::arg("grid_y"), py::arg("grid_z"),
             "Set grid dimensions")
        .def("set_order", &PMECoulomb::set_order,
             py::arg("order"),
             "Set B-spline interpolation order (4, 6, or 8 recommended)")
        .def_property_readonly("alpha", &PMECoulomb::alpha,
             "Get Ewald parameter")
        .def_property_readonly("order", &PMECoulomb::order,
             "Get B-spline order")
        .def_property_readonly("grid", &PMECoulomb::grid,
             "Get grid dimensions as (nx, ny, nz)")
        .def("set_charges", &PMECoulomb::set_charges,
             py::arg("charges"),
             "Set charges for all atoms (in units of e)")
        .def_property_readonly("charges", &PMECoulomb::charges,
             "Get charges")
        .def("cutoff", &PMECoulomb::cutoff,
             "Get real-space cutoff")
        .def("compute", &PMECoulomb::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial (requires 3D PBC)");

    // FMMCoulomb (Fast Multipole Method)
    py::class_<FMMCoulomb>(m, "FMMCoulomb")
        .def(py::init<int, int, int, int>(),
             py::arg("l_max") = 8, py::arg("n_level") = 3,
             py::arg("leaf_size") = 200, py::arg("periodic_images") = 1,
             "Construct FMM solver with expansion order, tree levels, max leaf size, and periodicity")
        .def("set_l_max", &FMMCoulomb::set_l_max,
             py::arg("l_max"),
             "Set maximum angular momentum for multipole expansion")
        .def("set_n_level", &FMMCoulomb::set_n_level,
             py::arg("n_level"),
             "Set number of tree levels")
        .def("set_leaf_size", &FMMCoulomb::set_leaf_size,
             py::arg("leaf_size"),
             "Set maximum particles per leaf")
        .def("set_periodic_images", &FMMCoulomb::set_periodic_images,
             py::arg("k"),
             "Set periodicity parameter (sum 3^k images)")
        .def_property_readonly("l_max", &FMMCoulomb::l_max,
             "Get maximum angular momentum")
        .def_property_readonly("n_level", &FMMCoulomb::n_level,
             "Get number of tree levels")
        .def_property_readonly("leaf_size", &FMMCoulomb::leaf_size,
             "Get maximum leaf size")
        .def("set_charges", &FMMCoulomb::set_charges,
             py::arg("charges"),
             "Set charges for all atoms (in units of e)")
        .def_property_readonly("charges", &FMMCoulomb::charges,
             "Get charges")
        .def("cutoff", &FMMCoulomb::cutoff,
             "Get cutoff (0 for FMM since it handles all interactions)")
        .def("compute", &FMMCoulomb::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial using FMM");
}
