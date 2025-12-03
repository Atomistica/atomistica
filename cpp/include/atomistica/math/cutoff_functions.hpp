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
#include <utility>

#include "../config.hpp"

namespace atomistica {

/**
 * @brief Result of cutoff function evaluation
 */
struct CutoffResult {
    Scalar fc;   // Cutoff function value
    Scalar dfc;  // Derivative with respect to r
};

/**
 * @brief Hard cutoff function (step function)
 *
 * fc(r) = 1 for r < cutoff, 0 otherwise
 */
class HardCutoff {
public:
    explicit HardCutoff(Scalar cutoff) : cutoff_(cutoff) {}

    Scalar cutoff() const { return cutoff_; }

    CutoffResult operator()(Scalar r) const {
        if (r < cutoff_) {
            return {1.0, 0.0};
        } else {
            return {0.0, 0.0};
        }
    }

private:
    Scalar cutoff_;
};

/**
 * @brief Shifted cutoff function
 *
 * fc(r) = 1 - r/cutoff for r < cutoff, 0 otherwise
 * Continuous but not smooth at cutoff
 */
class ShiftedCutoff {
public:
    explicit ShiftedCutoff(Scalar cutoff) : cutoff_(cutoff) {}

    Scalar cutoff() const { return cutoff_; }

    CutoffResult operator()(Scalar r) const {
        if (r < cutoff_) {
            return {1.0 - r / cutoff_, -1.0 / cutoff_};
        } else {
            return {0.0, 0.0};
        }
    }

private:
    Scalar cutoff_;
};

/**
 * @brief Trigonometric (cosine) cutoff function
 *
 * fc(r) = 0.5 * (1 + cos(pi * r / cutoff)) for r < cutoff
 * Smooth (continuous first derivative) at cutoff
 *
 * This is the standard cutoff used in Tersoff/Brenner potentials.
 */
class TrigonometricCutoff {
public:
    explicit TrigonometricCutoff(Scalar cutoff) : cutoff_(cutoff) {}

    Scalar cutoff() const { return cutoff_; }

    CutoffResult operator()(Scalar r) const {
        if (r >= cutoff_) {
            return {0.0, 0.0};
        }
        Scalar x = PI * r / cutoff_;
        Scalar fc = 0.5 * (1.0 + std::cos(x));
        Scalar dfc = -0.5 * PI / cutoff_ * std::sin(x);
        return {fc, dfc};
    }

private:
    Scalar cutoff_;
};

/**
 * @brief Inner/outer cutoff function for smooth transition region
 *
 * fc(r) = 1 for r < r_inner
 * fc(r) = 0.5 * (1 + cos(pi * (r - r_inner) / (r_outer - r_inner))) for r_inner <= r < r_outer
 * fc(r) = 0 for r >= r_outer
 */
class InnerOuterCutoff {
public:
    InnerOuterCutoff(Scalar inner, Scalar outer)
        : inner_(inner), outer_(outer), width_(outer - inner) {}

    Scalar cutoff() const { return outer_; }
    Scalar inner() const { return inner_; }

    CutoffResult operator()(Scalar r) const {
        if (r < inner_) {
            return {1.0, 0.0};
        } else if (r >= outer_) {
            return {0.0, 0.0};
        } else {
            Scalar x = PI * (r - inner_) / width_;
            Scalar fc = 0.5 * (1.0 + std::cos(x));
            Scalar dfc = -0.5 * PI / width_ * std::sin(x);
            return {fc, dfc};
        }
    }

private:
    Scalar inner_;
    Scalar outer_;
    Scalar width_;
};

/**
 * @brief Polynomial cutoff function
 *
 * fc(r) = (1 - (r/cutoff)^n)^m for r < cutoff
 * Default n=2, m=2 gives smooth behavior at both r=0 and r=cutoff
 */
template<int N = 2, int M = 2>
class PolynomialCutoff {
public:
    explicit PolynomialCutoff(Scalar cutoff) : cutoff_(cutoff) {}

    Scalar cutoff() const { return cutoff_; }

    CutoffResult operator()(Scalar r) const {
        if (r >= cutoff_) {
            return {0.0, 0.0};
        }
        Scalar x = r / cutoff_;
        Scalar xn = std::pow(x, N);
        Scalar term = 1.0 - xn;
        Scalar fc = std::pow(term, M);
        Scalar dfc = -M * std::pow(term, M - 1) * N * std::pow(x, N - 1) / cutoff_;
        return {fc, dfc};
    }

private:
    Scalar cutoff_;
};

// ============================================================================
// BOP-style cutoff functions (from Fortran atomistica)
// ============================================================================

/**
 * @brief Trigonometric "on" cutoff function
 *
 * Transition from 0 to 1:
 *   fc(r) = 0           for r <= r1
 *   fc(r) = 0.5*(1-cos(π*(r-r1)/(r2-r1)))  for r1 < r < r2
 *   fc(r) = 1           for r >= r2
 *
 * Differentiable once (C1).
 */
class TrigOnCutoff {
public:
    TrigOnCutoff() = default;

    TrigOnCutoff(Scalar r1, Scalar r2)
        : r1_(r1), r2_(r2), fac_(PI / (r2 - r1)) {}

    void init(Scalar r1, Scalar r2) {
        r1_ = r1;
        r2_ = r2;
        fac_ = PI / (r2_ - r1_);
    }

    Scalar r1() const { return r1_; }
    Scalar r2() const { return r2_; }

    CutoffResult operator()(Scalar r) const {
        if (r <= r1_) {
            return {0.0, 0.0};
        } else if (r >= r2_) {
            return {1.0, 0.0};
        } else {
            Scalar x = fac_ * (r - r1_);
            Scalar fc = 0.5 * (1.0 - std::cos(x));
            Scalar dfc = 0.5 * fac_ * std::sin(x);
            return {fc, dfc};
        }
    }

private:
    Scalar r1_ = 0.0;
    Scalar r2_ = 0.0;
    Scalar fac_ = 0.0;
};

/**
 * @brief Trigonometric "off" cutoff function
 *
 * Transition from 1 to 0:
 *   fc(r) = 1           for r <= r1
 *   fc(r) = 0.5*(1+cos(π*(r-r1)/(r2-r1)))  for r1 < r < r2
 *   fc(r) = 0           for r >= r2
 *
 * Differentiable once (C1). This is the standard BOP cutoff.
 */
class TrigOffCutoff {
public:
    TrigOffCutoff() = default;

    TrigOffCutoff(Scalar r1, Scalar r2)
        : r1_(r1), r2_(r2), fac_(PI / (r2 - r1)) {}

    void init(Scalar r1, Scalar r2) {
        r1_ = r1;
        r2_ = r2;
        fac_ = PI / (r2_ - r1_);
    }

    Scalar r1() const { return r1_; }
    Scalar r2() const { return r2_; }
    Scalar cutoff() const { return r2_; }

    CutoffResult operator()(Scalar r) const {
        if (r <= r1_) {
            return {1.0, 0.0};
        } else if (r >= r2_) {
            return {0.0, 0.0};
        } else {
            Scalar x = fac_ * (r - r1_);
            Scalar fc = 0.5 * (1.0 + std::cos(x));
            Scalar dfc = -0.5 * fac_ * std::sin(x);
            return {fc, dfc};
        }
    }

private:
    Scalar r1_ = 0.0;
    Scalar r2_ = 0.0;
    Scalar fac_ = 0.0;
};

/**
 * @brief Exponential cutoff function
 *
 * Based on f(x) = exp(-8*x^3), corrected so function, first and second
 * derivatives go to zero at x=1.
 *
 * Transition from 1 to 0:
 *   fc(r) = 1           for r <= r1
 *   fc(r) = corrected_exp(-8*x^3)  for r1 < r < r2, x = (r-r1)/(r2-r1)
 *   fc(r) = 0           for r >= r2
 *
 * Differentiable twice (C2). Used for screened potentials.
 */
class ExpCutoff {
public:
    ExpCutoff() = default;

    ExpCutoff(Scalar r1, Scalar r2) { init(r1, r2); }

    void init(Scalar r1, Scalar r2) {
        r1_ = r1;
        r2_ = r2;
        fac1_ = 1.0 / (r2_ - r1_);

        // Correction terms to ensure C2 continuity at x=1
        Scalar val1 = std::exp(-8.0);
        Scalar dval1 = -24.0 * val1;
        Scalar ddval1 = -48.0 * val1 - 24.0 * dval1;

        c_ = (-3.0 * dval1 + ddval1) / 3.0;
        d_ = (2.0 * dval1 - ddval1) / 4.0;
        fac2_ = 1.0 / (1.0 - val1 - c_ - d_);
        off_ = val1 + c_ + d_;
    }

    Scalar r1() const { return r1_; }
    Scalar r2() const { return r2_; }
    Scalar cutoff() const { return r2_; }

    CutoffResult operator()(Scalar r) const {
        if (r <= r1_) {
            return {1.0, 0.0};
        } else if (r >= r2_) {
            return {0.0, 0.0};
        } else {
            Scalar x = fac1_ * (r - r1_);
            Scalar x2 = x * x;
            Scalar x3 = x * x2;

            Scalar exp_val = std::exp(-8.0 * x3);
            Scalar dexp_val = -24.0 * x2 * exp_val;

            // Apply correction for C2 continuity
            Scalar fc = fac2_ * (exp_val + c_ * x3 + d_ * x2 * x2 - off_);
            Scalar dfc = fac1_ * fac2_ * (dexp_val + 3.0 * c_ * x2 + 4.0 * d_ * x3);

            return {fc, dfc};
        }
    }

private:
    Scalar r1_ = 0.0;
    Scalar r2_ = 0.0;
    Scalar fac1_ = 0.0;
    Scalar fac2_ = 0.0;
    Scalar c_ = 0.0;
    Scalar d_ = 0.0;
    Scalar off_ = 0.0;
};

/**
 * @brief Polymorphic cutoff function wrapper
 *
 * Can hold either TrigOffCutoff or ExpCutoff for runtime selection.
 */
class BOPCutoff {
public:
    enum class Type { TrigOff, Exp };

    BOPCutoff() = default;

    BOPCutoff(Type type, Scalar r1, Scalar r2) : type_(type) {
        if (type == Type::TrigOff) {
            trig_off_.init(r1, r2);
        } else {
            exp_.init(r1, r2);
        }
    }

    void init(Type type, Scalar r1, Scalar r2) {
        type_ = type;
        if (type == Type::TrigOff) {
            trig_off_.init(r1, r2);
        } else {
            exp_.init(r1, r2);
        }
    }

    Scalar r1() const {
        return (type_ == Type::TrigOff) ? trig_off_.r1() : exp_.r1();
    }

    Scalar r2() const {
        return (type_ == Type::TrigOff) ? trig_off_.r2() : exp_.r2();
    }

    Scalar cutoff() const { return r2(); }

    CutoffResult operator()(Scalar r) const {
        if (type_ == Type::TrigOff) {
            return trig_off_(r);
        } else {
            return exp_(r);
        }
    }

private:
    Type type_ = Type::TrigOff;
    TrigOffCutoff trig_off_;
    ExpCutoff exp_;
};

} // namespace atomistica
