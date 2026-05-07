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
#include <atomistica/potentials/bop/juslin.hpp>
#include <atomistica/potentials/pair/simple_pairs.hpp>
#include <atomistica/tightbinding/tightbinding.hpp>

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
        return std::vector<std::string>{
            "Tersoff_PRB_39_5566_Si_C",
            "Goumri_Said_ChemPhys_302_135_Al_N",
            "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N",
            "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr",
        };
    }, "List available built-in Tersoff parameter sets");

    // =========================================================================
    // Brenner Potential
    // =========================================================================

    // BrennerElementParams (needed before Brenner class for default argument)
    py::class_<BrennerElementParams>(m, "BrennerElementParams")
        .def(py::init<>());

    // BrennerPairParams
    py::class_<BrennerPairParams>(m, "BrennerPairParams")
        .def(py::init<>())
        .def_readwrite("D0", &BrennerPairParams::D0)
        .def_readwrite("r0", &BrennerPairParams::r0)
        .def_readwrite("S", &BrennerPairParams::S)
        .def_readwrite("beta", &BrennerPairParams::beta)
        .def_readwrite("gamma", &BrennerPairParams::gamma)
        .def_readwrite("c", &BrennerPairParams::c)
        .def_readwrite("d", &BrennerPairParams::d)
        .def_readwrite("h", &BrennerPairParams::h)
        .def_readwrite("mu", &BrennerPairParams::mu)
        .def_readwrite("n", &BrennerPairParams::n)
        .def_readwrite("m", &BrennerPairParams::m)
        .def_readwrite("r1", &BrennerPairParams::r1)
        .def_readwrite("r2", &BrennerPairParams::r2)
        .def("precompute", &BrennerPairParams::precompute);

    // Brenner potential (non-screened)
    py::class_<Brenner<false>>(m, "Brenner")
        .def(py::init<>())
        .def("add_element", &Brenner<false>::add_element,
             py::arg("Z"), py::arg("params") = BrennerElementParams{},
             "Add element with given atomic number")
        .def("set_pair_params", &Brenner<false>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"),
             "Set pair parameters for element pair")
        .def("load_parameters", &Brenner<false>::load_parameters,
             py::arg("name"),
             "Load built-in parameter set by name")
        .def("cutoff", &Brenner<false>::cutoff,
             "Get maximum cutoff radius")
        .def("num_elements", &Brenner<false>::num_elements,
             "Get number of elements defined")
        .def("element_index", &Brenner<false>::element_index,
             py::arg("Z"),
             "Get internal element index for atomic number Z (-1 if not found)")
        .def("pair_type", &Brenner<false>::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type index for element pair")
        .def("compute", &Brenner<false>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        // Expose internal functions for testing/analysis
        .def("repulsive", [](const Brenner<false>& pot, int ptype, Scalar r) {
            auto [V, dV] = pot.repulsive(ptype, r);
            return std::make_pair(V, dV);
        }, py::arg("ptype"), py::arg("r"))
        .def("attractive", [](const Brenner<false>& pot, int ptype, Scalar r) {
            auto [V, dV] = pot.attractive(ptype, r);
            return std::make_pair(V, dV);
        }, py::arg("ptype"), py::arg("r"))
        .def("angular_function", [](const Brenner<false>& pot,
                int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) {
            auto [g, dg] = pot.angular_function(eli, elj, elk, ptype_ij, ptype_ik, cos_theta);
            return std::make_pair(g, dg);
        }, py::arg("eli"), py::arg("elj"), py::arg("elk"),
           py::arg("ptype_ij"), py::arg("ptype_ik"), py::arg("cos_theta"))
        .def("bond_order", [](const Brenner<false>& pot, int eli, int ptype, Scalar z) {
            auto [b, db] = pot.bond_order(eli, ptype, z);
            return std::make_pair(b, db);
        }, py::arg("eli"), py::arg("ptype"), py::arg("z"));

    // Screened Brenner potential
    py::class_<Brenner<true>>(m, "BrennerScr")
        .def(py::init<>())
        .def("add_element", &Brenner<true>::add_element,
             py::arg("Z"), py::arg("params") = BrennerElementParams{},
             "Add element with given atomic number")
        .def("set_pair_params", &Brenner<true>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"),
             "Set pair parameters for element pair")
        .def("load_parameters", &Brenner<true>::load_parameters,
             py::arg("name"),
             "Load built-in parameter set by name")
        .def("cutoff", &Brenner<true>::cutoff,
             "Get maximum cutoff radius")
        .def("num_elements", &Brenner<true>::num_elements,
             "Get number of elements defined")
        .def("element_index", &Brenner<true>::element_index,
             py::arg("Z"),
             "Get internal element index for atomic number Z (-1 if not found)")
        .def("pair_type", &Brenner<true>::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type index for element pair")
        .def("compute", &Brenner<true>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        .def("screening_params", &Brenner<true>::screening_params,
             py::arg("ptype"),
             "Get screening parameters for pair type",
             py::return_value_policy::reference_internal);

    // Available Brenner parameter sets
    m.def("available_brenner_parameters", []() {
        return std::vector<std::string>{
            "Erhart_PRB_71_035211_SiC",
            "Albe_PRB_65_195124_PtC",
            "Henriksson_PRB_79_144107_FeC",
            "Kioseoglou_PSSb_245_1118_AlN",
            "Brenner_PRB_42_9458_C_I",
            "Brenner_PRB_42_9458_C_II",
        };
    }, "List available built-in Brenner parameter sets");

    // =========================================================================
    // Kumagai Potential
    // =========================================================================

    // KumagaiElementParams (needed before Kumagai class for default argument)
    py::class_<KumagaiElementParams>(m, "KumagaiElementParams")
        .def(py::init<>())
        .def_readwrite("c1", &KumagaiElementParams::c1)
        .def_readwrite("c2", &KumagaiElementParams::c2)
        .def_readwrite("c3", &KumagaiElementParams::c3)
        .def_readwrite("c4", &KumagaiElementParams::c4)
        .def_readwrite("c5", &KumagaiElementParams::c5)
        .def_readwrite("h", &KumagaiElementParams::h)
        .def_readwrite("eta", &KumagaiElementParams::eta)
        .def_readwrite("delta", &KumagaiElementParams::delta)
        .def("precompute", &KumagaiElementParams::precompute);

    // KumagaiPairParams
    py::class_<KumagaiPairParams>(m, "KumagaiPairParams")
        .def(py::init<>())
        .def_readwrite("A", &KumagaiPairParams::A)
        .def_readwrite("B", &KumagaiPairParams::B)
        .def_readwrite("lambda_", &KumagaiPairParams::lambda)
        .def_readwrite("mu", &KumagaiPairParams::mu)
        .def_readwrite("alpha", &KumagaiPairParams::alpha)
        .def_readwrite("beta", &KumagaiPairParams::beta)
        .def_readwrite("r1", &KumagaiPairParams::r1)
        .def_readwrite("r2", &KumagaiPairParams::r2)
        .def("precompute", &KumagaiPairParams::precompute);

    // Kumagai potential (non-screened)
    py::class_<Kumagai<false>>(m, "Kumagai")
        .def(py::init<>())
        .def("add_element", &Kumagai<false>::add_element,
             py::arg("Z"), py::arg("params") = KumagaiElementParams{},
             "Add element with given atomic number")
        .def("set_pair_params", &Kumagai<false>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"),
             "Set pair parameters for element pair")
        .def("load_parameters", &Kumagai<false>::load_parameters,
             py::arg("name"),
             "Load built-in parameter set by name")
        .def("cutoff", &Kumagai<false>::cutoff,
             "Get maximum cutoff radius")
        .def("num_elements", &Kumagai<false>::num_elements,
             "Get number of elements defined")
        .def("element_index", &Kumagai<false>::element_index,
             py::arg("Z"),
             "Get internal element index for atomic number Z (-1 if not found)")
        .def("pair_type", &Kumagai<false>::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type index for element pair")
        .def("compute", &Kumagai<false>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        // Expose internal functions for testing/analysis
        .def("repulsive", [](const Kumagai<false>& pot, int ptype, Scalar r) {
            auto [V, dV] = pot.repulsive(ptype, r);
            return std::make_pair(V, dV);
        }, py::arg("ptype"), py::arg("r"))
        .def("attractive", [](const Kumagai<false>& pot, int ptype, Scalar r) {
            auto [V, dV] = pot.attractive(ptype, r);
            return std::make_pair(V, dV);
        }, py::arg("ptype"), py::arg("r"))
        .def("angular_function", [](const Kumagai<false>& pot,
                int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) {
            auto [g, dg] = pot.angular_function(eli, elj, elk, ptype_ij, ptype_ik, cos_theta);
            return std::make_pair(g, dg);
        }, py::arg("eli"), py::arg("elj"), py::arg("elk"),
           py::arg("ptype_ij"), py::arg("ptype_ik"), py::arg("cos_theta"))
        .def("bond_order", [](const Kumagai<false>& pot, int eli, int ptype, Scalar z) {
            auto [b, db] = pot.bond_order(eli, ptype, z);
            return std::make_pair(b, db);
        }, py::arg("eli"), py::arg("ptype"), py::arg("z"));

    // Screened Kumagai potential
    py::class_<Kumagai<true>>(m, "KumagaiScr")
        .def(py::init<>())
        .def("add_element", &Kumagai<true>::add_element,
             py::arg("Z"), py::arg("params") = KumagaiElementParams{},
             "Add element with given atomic number")
        .def("set_pair_params", &Kumagai<true>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"),
             "Set pair parameters for element pair")
        .def("load_parameters", &Kumagai<true>::load_parameters,
             py::arg("name"),
             "Load built-in parameter set by name")
        .def("cutoff", &Kumagai<true>::cutoff,
             "Get maximum cutoff radius")
        .def("num_elements", &Kumagai<true>::num_elements,
             "Get number of elements defined")
        .def("element_index", &Kumagai<true>::element_index,
             py::arg("Z"),
             "Get internal element index for atomic number Z (-1 if not found)")
        .def("pair_type", &Kumagai<true>::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type index for element pair")
        .def("compute", &Kumagai<true>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        .def("screening_params", &Kumagai<true>::screening_params,
             py::arg("ptype"),
             "Get screening parameters for pair type",
             py::return_value_policy::reference_internal);

    // Available Kumagai parameter sets
    m.def("available_kumagai_parameters", []() {
        return std::vector<std::string>{
            "Kumagai_CompMaterSci_39_457_Si"
        };
    }, "List available built-in Kumagai parameter sets");

    // =========================================================================
    // REBO2 (2nd Generation Brenner Potential)
    // =========================================================================

    py::class_<REBO2>(m, "REBO2")
        .def(py::init<>())
        .def("load_default_parameters", &REBO2::load_default_parameters,
             "Load default Brenner 2002 parameters")
        .def("cutoff", &REBO2::cutoff,
             "Get cutoff radius")
        .def("element_type", &REBO2::element_type,
             py::arg("Z"),
             "Get internal element type from atomic number (1=C, 3=H, or -1 if unsupported)")
        .def_static("pair_type", &REBO2::pair_type,
             py::arg("eli"), py::arg("elj"),
             "Get pair type from two element types")
        .def("repulsive", [](const REBO2& pot, int ptype, Scalar r) {
            auto [val, deriv] = pot.repulsive(ptype, r);
            return std::make_pair(val, deriv);
        }, py::arg("ptype"), py::arg("r"),
           "Evaluate repulsive pair function V_R(r)")
        .def("attractive", [](const REBO2& pot, int ptype, Scalar r) {
            auto [val, deriv] = pot.attractive(ptype, r);
            return std::make_pair(val, deriv);
        }, py::arg("ptype"), py::arg("r"),
           "Evaluate attractive pair function V_A(r)")
        .def("angular_function", [](const REBO2& pot, int el_type, Scalar cos_theta, Scalar N) {
            auto [g, dg_dcos, dg_dN] = pot.angular_function(el_type, cos_theta, N);
            return std::make_tuple(g, dg_dcos, dg_dN);
        }, py::arg("el_type"), py::arg("cos_theta"), py::arg("N"),
           "Evaluate angular function g(cos_theta, N)")
        .def("bond_order_func", [](const REBO2& pot, int el_type, Scalar z) {
            auto [b, db] = pot.bond_order_func(el_type, z);
            return std::make_pair(b, db);
        }, py::arg("el_type"), py::arg("z"),
           "Evaluate bond order function b(1+z)")
        .def("compute", &REBO2::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute energy, forces, and virial")
        // Expose parameters for modification/inspection
        .def_readwrite("cc_B1", &REBO2::cc_B1)
        .def_readwrite("cc_B2", &REBO2::cc_B2)
        .def_readwrite("cc_B3", &REBO2::cc_B3)
        .def_readwrite("cc_beta1", &REBO2::cc_beta1)
        .def_readwrite("cc_beta2", &REBO2::cc_beta2)
        .def_readwrite("cc_beta3", &REBO2::cc_beta3)
        .def_readwrite("cc_Q", &REBO2::cc_Q)
        .def_readwrite("cc_A", &REBO2::cc_A)
        .def_readwrite("cc_alpha", &REBO2::cc_alpha)
        .def_readwrite("ch_B1", &REBO2::ch_B1)
        .def_readwrite("ch_beta1", &REBO2::ch_beta1)
        .def_readwrite("ch_Q", &REBO2::ch_Q)
        .def_readwrite("ch_A", &REBO2::ch_A)
        .def_readwrite("ch_alpha", &REBO2::ch_alpha)
        .def_readwrite("hh_B1", &REBO2::hh_B1)
        .def_readwrite("hh_beta1", &REBO2::hh_beta1)
        .def_readwrite("hh_Q", &REBO2::hh_Q)
        .def_readwrite("hh_A", &REBO2::hh_A)
        .def_readwrite("hh_alpha", &REBO2::hh_alpha);

    // REBO2 pair type constants
    m.attr("REBO2_C_C") = REBO2_C_C;
    m.attr("REBO2_C_H") = REBO2_C_H;
    m.attr("REBO2_H_H") = REBO2_H_H;
    m.attr("REBO2_C") = REBO2_C;
    m.attr("REBO2_H") = REBO2_H;

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

    // ========================================================================
    // Tight-Binding / DFTB
    // ========================================================================

    // TBElementParams
    py::class_<tb::TBElementParams>(m, "TBElementParams")
        .def(py::init<>())
        .def_readwrite("symbol", &tb::TBElementParams::symbol)
        .def_readwrite("atomic_number", &tb::TBElementParams::atomic_number)
        .def_readwrite("num_orbitals", &tb::TBElementParams::num_orbitals)
        .def_readwrite("l_max", &tb::TBElementParams::l_max)
        .def_readwrite("hubbard_U", &tb::TBElementParams::hubbard_U)
        .def_readwrite("valence_electrons", &tb::TBElementParams::valence_electrons)
        .def_property("onsite",
            [](const tb::TBElementParams& e) {
                std::vector<Scalar> v(e.onsite.begin(), e.onsite.begin() + e.num_orbitals);
                return v;
            },
            [](tb::TBElementParams& e, const std::vector<Scalar>& v) {
                for (size_t i = 0; i < v.size() && i < 9; ++i) {
                    e.onsite[i] = v[i];
                }
            })
        .def("is_s_only", &tb::TBElementParams::is_s_only)
        .def("is_sp", &tb::TBElementParams::is_sp)
        .def("is_spd", &tb::TBElementParams::is_spd);

    // Predefined element parameters
    m.def("carbon_mio", &tb::parameters::carbon_mio, "Carbon parameters for mio-1-1 DFTB");
    m.def("hydrogen_mio", &tb::parameters::hydrogen_mio, "Hydrogen parameters for mio-1-1 DFTB");
    m.def("oxygen_mio", &tb::parameters::oxygen_mio, "Oxygen parameters for mio-1-1 DFTB");
    m.def("nitrogen_mio", &tb::parameters::nitrogen_mio, "Nitrogen parameters for mio-1-1 DFTB");

    // SCCParams
    py::class_<tb::SCCParams>(m, "SCCParams")
        .def(py::init<>())
        .def_readwrite("max_iterations", &tb::SCCParams::max_iterations)
        .def_readwrite("convergence_threshold", &tb::SCCParams::convergence_threshold)
        .def_readwrite("mixing_parameter", &tb::SCCParams::mixing_parameter)
        .def_readwrite("anderson_memory", &tb::SCCParams::anderson_memory)
        .def_readwrite("enable_dftb3", &tb::SCCParams::enable_dftb3)
        .def_readwrite("zeta", &tb::SCCParams::zeta);

    // SolverParams
    py::class_<tb::SolverParams>(m, "SolverParams")
        .def(py::init<>())
        .def_readwrite("electronic_temperature", &tb::SolverParams::electronic_temperature)
        .def_readwrite("use_divide_and_conquer", &tb::SolverParams::use_divide_and_conquer);

    // DenseHamiltonian - for accessing internal data
    py::class_<tb::DenseHamiltonian>(m, "DenseHamiltonian")
        .def_readonly("num_atoms", &tb::DenseHamiltonian::num_atoms)
        .def_readonly("num_orbitals", &tb::DenseHamiltonian::num_orbitals)
        .def_property_readonly("H", [](const tb::DenseHamiltonian& h) { return h.H; })
        .def_property_readonly("S", [](const tb::DenseHamiltonian& h) { return h.S; })
        .def_property_readonly("rho", [](const tb::DenseHamiltonian& h) { return h.rho; })
        .def_property_readonly("eigenvalues", [](const tb::DenseHamiltonian& h) { return h.eigenvalues; })
        .def_property_readonly("eigenvectors", [](const tb::DenseHamiltonian& h) { return h.eigenvectors; })
        .def_property_readonly("occupation", [](const tb::DenseHamiltonian& h) { return h.occupation; })
        .def_property_readonly("charges", [](const tb::DenseHamiltonian& h) { return h.charges; })
        .def_readonly("band_energy", &tb::DenseHamiltonian::band_energy)
        .def_readonly("repulsive_energy", &tb::DenseHamiltonian::repulsive_energy)
        .def_readonly("fermi_level", &tb::DenseHamiltonian::fermi_level);

    // MaterialsDatabase
    py::class_<tb::MaterialsDatabase>(m, "MaterialsDatabase")
        .def(py::init<>())
        .def("load_skf_directory", &tb::MaterialsDatabase::load_skf_directory,
             py::arg("path"),
             "Load SKF files from directory")
        .def("add_element", &tb::MaterialsDatabase::add_element,
             py::arg("elem"),
             "Add element parameters manually")
        .def("has_element", &tb::MaterialsDatabase::has_element,
             py::arg("Z"),
             "Check if element exists in database")
        .def("get_element", &tb::MaterialsDatabase::get_element,
             py::arg("Z"),
             "Get element parameters by atomic number")
        .def("load_pair", &tb::MaterialsDatabase::load_pair,
             py::arg("Z1"), py::arg("Z2"),
             "Load pair parameters from SKF file")
        .def("get_cutoff", &tb::MaterialsDatabase::get_cutoff,
             py::arg("Z1"), py::arg("Z2"),
             "Get cutoff for pair")
        .def("get_max_cutoff", &tb::MaterialsDatabase::get_max_cutoff,
             "Get maximum cutoff across all loaded pairs");

    // DFTB potential
    py::class_<tb::DFTB>(m, "DFTB")
        .def(py::init<const std::string&, bool>(),
             py::arg("skf_path") = "",
             py::arg("enable_scc") = false,
             "Create DFTB potential with optional SKF path and SCC")
        .def("name", &tb::DFTB::name,
             "Get potential name")
        .def("cutoff", &tb::DFTB::cutoff,
             "Get cutoff distance")
        .def("add_element", &tb::DFTB::add_element,
             py::arg("elem"),
             "Add element to materials database")
        .def("load_pair", &tb::DFTB::load_pair,
             py::arg("Z1"), py::arg("Z2"),
             "Load pair parameters from SKF file")
        .def("set_skf_path", &tb::DFTB::set_skf_path,
             py::arg("path"),
             "Set SKF directory path")
        .def("set_scc", &tb::DFTB::set_scc,
             py::arg("enable"),
             "Enable/disable SCC")
        .def("set_scc_params", &tb::DFTB::set_scc_params,
             py::arg("params"),
             "Set SCC parameters")
        .def("set_solver_params", &tb::DFTB::set_solver_params,
             py::arg("params"),
             "Set solver parameters")
        .def("init", &tb::DFTB::init,
             py::arg("system"),
             "Initialize potential for atomic system")
        .def("compute_energy", &tb::DFTB::compute_energy,
             py::arg("system"), py::arg("neighbors"),
             "Compute energy only")
        .def("compute", [](tb::DFTB& dftb, AtomicSystem& system, NeighborList& neighbors,
                          bool compute_forces, bool compute_virial) {
                 MatX3 forces;
                 Scalar energy = dftb.compute(system, neighbors, forces);

                 PotentialResults results;
                 results.energy = energy;

                 if (compute_forces) {
                     // Copy forces to system
                     for (std::size_t i = 0; i < static_cast<std::size_t>(system.num_atoms()); ++i) {
                         system.forces().col(i) = forces.row(i).transpose();
                     }
                 }

                 return results;
             },
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = false,
             "Compute energy and forces")
        .def_property_readonly("hamiltonian",
             [](tb::DFTB& dftb) -> const tb::DenseHamiltonian& { return dftb.hamiltonian(); },
             py::return_value_policy::reference_internal,
             "Get Hamiltonian data structure")
        .def_property_readonly("eigenvalues", &tb::DFTB::eigenvalues,
             "Get eigenvalues")
        .def_property_readonly("fermi_level", &tb::DFTB::fermi_level,
             "Get Fermi level")
        .def_property_readonly("charges", &tb::DFTB::charges,
             "Get Mulliken charges")
        .def_property_readonly("band_energy", &tb::DFTB::band_energy,
             "Get band energy")
        .def_property_readonly("repulsive_energy", &tb::DFTB::repulsive_energy,
             "Get repulsive energy")
        .def_property_readonly("materials",
             [](tb::DFTB& dftb) -> tb::MaterialsDatabase& { return dftb.materials(); },
             py::return_value_policy::reference_internal,
             "Get materials database");

    // Slater-Koster transformation function (for testing/debugging)
    m.def("sk_transform", &tb::transform_orb,
          py::arg("a"), py::arg("b"), py::arg("c"), py::arg("sk"),
          "Transform SK integrals to Cartesian matrix element");

    m.def("sk_transform_derivative", &tb::transform_orb_derivative,
          py::arg("a"), py::arg("b"), py::arg("c"), py::arg("r"),
          py::arg("sk"), py::arg("dsk"),
          "Compute derivative of SK-transformed matrix element");

    // =========================================================================
    // Juslin Potential
    // =========================================================================

    // Juslin (non-screened)
    py::class_<Juslin<false>>(m, "Juslin")
        .def(py::init<>())
        .def("add_element", &Juslin<false>::add_element,
             py::arg("Z"), py::arg("params") = BrennerElementParams{})
        .def("set_pair_params", &Juslin<false>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"))
        .def("set_triplet_params", &Juslin<false>::set_triplet_params,
             py::arg("eli"), py::arg("elj"), py::arg("elk"),
             py::arg("alpha"), py::arg("omega"), py::arg("m"))
        .def("load_parameters", &Juslin<false>::load_parameters,
             py::arg("name"), "Load built-in parameter set by name")
        .def("cutoff", &Juslin<false>::cutoff)
        .def("num_elements", &Juslin<false>::num_elements)
        .def("element_index", &Juslin<false>::element_index, py::arg("Z"))
        .def("pair_type", &Juslin<false>::pair_type, py::arg("eli"), py::arg("elj"))
        .def("compute", &Juslin<false>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    // Screened Juslin
    py::class_<Juslin<true>>(m, "JuslinScr")
        .def(py::init<>())
        .def("add_element", &Juslin<true>::add_element,
             py::arg("Z"), py::arg("params") = BrennerElementParams{})
        .def("set_pair_params", &Juslin<true>::set_pair_params,
             py::arg("Z1"), py::arg("Z2"), py::arg("params"))
        .def("set_triplet_params", &Juslin<true>::set_triplet_params,
             py::arg("eli"), py::arg("elj"), py::arg("elk"),
             py::arg("alpha"), py::arg("omega"), py::arg("m"))
        .def("load_parameters", &Juslin<true>::load_parameters,
             py::arg("name"), "Load built-in parameter set by name")
        .def("cutoff", &Juslin<true>::cutoff)
        .def("num_elements", &Juslin<true>::num_elements)
        .def("element_index", &Juslin<true>::element_index, py::arg("Z"))
        .def("pair_type", &Juslin<true>::pair_type, py::arg("eli"), py::arg("elj"))
        .def("compute", &Juslin<true>::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    m.def("available_juslin_parameters", []() {
        return std::vector<std::string>{"Juslin_JAP_98_123520_WCH"};
    });

    // =========================================================================
    // REBO2Scr
    // =========================================================================

    py::class_<REBO2Scr>(m, "REBO2Scr")
        .def(py::init<>())
        .def("load_default_parameters", &REBO2Scr::load_default_parameters)
        .def("cutoff", &REBO2Scr::cutoff)
        .def("element_type", &REBO2Scr::element_type, py::arg("Z"))
        .def_static("pair_type", &REBO2Scr::pair_type, py::arg("eli"), py::arg("elj"))
        .def("repulsive", [](const REBO2Scr& pot, int ptype, Scalar r) {
            auto [val, deriv] = pot.repulsive(ptype, r);
            return std::make_pair(val, deriv);
        }, py::arg("ptype"), py::arg("r"))
        .def("attractive", [](const REBO2Scr& pot, int ptype, Scalar r) {
            auto [val, deriv] = pot.attractive(ptype, r);
            return std::make_pair(val, deriv);
        }, py::arg("ptype"), py::arg("r"))
        .def("compute", &REBO2Scr::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    // =========================================================================
    // Simple pair potentials
    // =========================================================================

    py::class_<BornMayer>(m, "BornMayer")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar, Scalar, int, int>(),
             py::arg("A"), py::arg("rho"), py::arg("cutoff"),
             py::arg("Z1") = 0, py::arg("Z2") = 0)
        .def_readwrite("A", &BornMayer::A)
        .def_readwrite("rho", &BornMayer::rho)
        .def_readwrite("cutoff_radius", &BornMayer::cutoff_radius)
        .def_readwrite("Z1", &BornMayer::Z1)
        .def_readwrite("Z2", &BornMayer::Z2)
        .def("cutoff", &BornMayer::cutoff)
        .def("bind_to", &BornMayer::bind_to)
        .def("compute", &BornMayer::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    py::class_<Harmonic>(m, "Harmonic")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar, Scalar, bool, int, int>(),
             py::arg("k"), py::arg("r0"), py::arg("cutoff"),
             py::arg("shift") = false, py::arg("Z1") = 0, py::arg("Z2") = 0)
        .def_readwrite("k", &Harmonic::k)
        .def_readwrite("r0", &Harmonic::r0)
        .def_readwrite("cutoff_radius", &Harmonic::cutoff_radius)
        .def_readwrite("shift", &Harmonic::shift)
        .def_readwrite("Z1", &Harmonic::Z1)
        .def_readwrite("Z2", &Harmonic::Z2)
        .def("cutoff", &Harmonic::cutoff)
        .def("bind_to", &Harmonic::bind_to)
        .def("compute", &Harmonic::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    py::class_<DoubleHarmonic>(m, "DoubleHarmonic")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar, Scalar, Scalar, Scalar, int, int>(),
             py::arg("k1"), py::arg("r1"), py::arg("k2"), py::arg("r2"),
             py::arg("cutoff"), py::arg("Z1") = 0, py::arg("Z2") = 0)
        .def_readwrite("k1", &DoubleHarmonic::k1)
        .def_readwrite("r1", &DoubleHarmonic::r1)
        .def_readwrite("k2", &DoubleHarmonic::k2)
        .def_readwrite("r2", &DoubleHarmonic::r2)
        .def_readwrite("cutoff_radius", &DoubleHarmonic::cutoff_radius)
        .def_readwrite("Z1", &DoubleHarmonic::Z1)
        .def_readwrite("Z2", &DoubleHarmonic::Z2)
        .def("cutoff", &DoubleHarmonic::cutoff)
        .def("bind_to", &DoubleHarmonic::bind_to)
        .def("compute", &DoubleHarmonic::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    py::class_<R6>(m, "R6")
        .def(py::init<>())
        .def(py::init<Scalar, Scalar, Scalar, int, int>(),
             py::arg("A"), py::arg("r0"), py::arg("cutoff"),
             py::arg("Z1") = 0, py::arg("Z2") = 0)
        .def_readwrite("A", &R6::A)
        .def_readwrite("r0", &R6::r0)
        .def_readwrite("cutoff_radius", &R6::cutoff_radius)
        .def_readwrite("Z1", &R6::Z1)
        .def_readwrite("Z2", &R6::Z2)
        .def("cutoff", &R6::cutoff)
        .def("bind_to", &R6::bind_to)
        .def("compute", &R6::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true, py::arg("compute_virial") = true);

    // =========================================================================
    // DFT-D3 dispersion correction
    // =========================================================================

    py::class_<DFTD3Disp>(m, "DFTD3Disp")
        .def(py::init<>())
        .def_readwrite("s6",           &DFTD3Disp::s6)
        .def_readwrite("s8",           &DFTD3Disp::s8)
        .def_readwrite("a1",           &DFTD3Disp::a1)
        .def_readwrite("a2",           &DFTD3Disp::a2)
        .def_readwrite("sr6",          &DFTD3Disp::sr6)
        .def_readwrite("sr8",          &DFTD3Disp::sr8)
        .def_readwrite("alpha6",       &DFTD3Disp::alpha6)
        .def_readwrite("cutoff_radius",&DFTD3Disp::cutoff_radius)
        .def_readwrite("cutoff_cn",    &DFTD3Disp::cutoff_cn)
        .def_readwrite("use_bj_damping",&DFTD3Disp::use_bj_damping)
        .def("cutoff",  &DFTD3Disp::cutoff_impl,
             "Interaction cutoff radius (Å)")
        .def("bind_to", &DFTD3Disp::bind_to,
             py::arg("system"), py::arg("neighbors"),
             "Initialise the DFT-D3 calculator for the given system")
        .def("compute", &DFTD3Disp::compute,
             py::arg("system"), py::arg("neighbors"),
             py::arg("compute_forces") = true,
             py::arg("compute_virial") = true,
             "Compute DFT-D3 dispersion energy, forces and virial")
        .def_property_readonly("have_sdftd3", [](const DFTD3Disp&) {
#ifdef HAVE_SDFTD3
            return true;
#else
            return false;
#endif
        }, "True if the s-dftd3 library was found at build time");
}
