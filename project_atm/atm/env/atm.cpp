#include "atm.h"
#include <cmath>
#include <stdexcept>

/*
 * Atmosphere model intended for aircraft trajectory prediction.
 * Provides atmospheric properties based on geopotential altitude based on
 * temperature and pressure offsets.
 * Coincides with ICAO Standard Atmosphere (ISA) when both offsets are zero.
 * Copyright (C) <2020r>  <Eduardo Gallo>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Author contact information: edugallo@yahoo.com
*/

// CLASS ATM
// =========
// =========

/* ===== ===== ===== Atmospheric constants ===== ===== ===== */
const double	env::atm::_R = 287.05287;
/* air specific constant [J/(kg*degK)] = [m2/(degK*s2)] */
const double	env::atm::_kappa = 1.4;
/* air adiabatic index [-] */
const double	env::atm::_betaT = -6.5e-3;
/* troposphere thermal gradient [degK/m] */
const double	env::atm::_g0_mps2 = 9.80665;
/* standard acceleration of free fall [m/s2] */
const double    env::atm::_g0_R = - env::atm::_g0_mps2 / env::atm::_R;
/* division of g0 by R, with negative sign [degK/m] */
const double    env::atm::_Hp1_m(11000.);
/* tropopause geopotential altitude */

/* ===== ===== ===== Standard mean sea level conditions ===== ===== ==== */
/* Conditions at standard mean sea level Hp = 0 if both temperature and pressure offsets are zero */
const double	env::atm::_p0_pa(101325.);
/* standard mean sea level pressure */
const double	env::atm::_T0_degK(288.15);
/* standard mean sea level temperature */

/* ===== ===== ===== Constructors ===== ===== ===== */
env::atm::atm(const double& DeltaT_degK, const double& Deltap_pa)
: _DeltaT_degK(DeltaT_degK), _Deltap_pa(Deltap_pa) {
}
/* constructor based on temperature and pressure differentials at mean sea level. */

env::atm::atm(const double& T_degK_datum, const double& p_pa_datum, const double& H_m_datum)
: _DeltaT_degK(0.), _Deltap_pa(0.) {
    const double Hp_m_datum = env::atm::p2Hp(p_pa_datum); // pressure altitude
    if (Hp_m_datum > _Hp1_m) { // Hpd shall be in troposphere
        throw std::runtime_error("Input pressure shall be below tropopause");
    }
    const double Tisa_degK_datum = env::atm::Hp2Tisa(Hp_m_datum); // standard temperature
    const double DeltaT_degK = T_degK_datum - Tisa_degK_datum; // temperature offset

    constructor_aux Func(this);
    this->set(DeltaT_degK, 0.); // object with proper _DeltaT created, used by acm_aux functions
    std::vector<double> par(2);
    par[0] = Tisa_degK_datum;
    par[1] = H_m_datum;

    double Tisa_degK_msl = Func.find_zero_secant(par, 0, _T0_degK - 30.0, _T0_degK + 30.0, 1e-13, 100); // standard temperature at mean sea level
    double Hp_m_msl	= env::atm::Tisa2Hp(Tisa_degK_msl); // pressure altitude at mean sea level
    double p_pa_msl = env::atm::Hp2p(Hp_m_msl); // pressure at mean sea level
    double Deltap_pa = p_pa_msl - _p0_pa; // pressure offset
    this->set(DeltaT_degK, Deltap_pa); // fill up temperature and pressure offsets
}
/* constructor based on given point temperature, pressure, and geopotential altitude. This point, normally at
 * airport or meteorological station, is assumed to be in the troposphere). */

void env::atm::set(const double& DeltaT_degK, const double& Deltap_pa) {
    _DeltaT_degK = DeltaT_degK;
    _Deltap_pa = Deltap_pa;
}
/* sets the temperature and pressure offsets at mean sea level */

/* ===== ===== ===== Static Functions ===== ===== ===== */
/* ==================================================== */
double env::atm::Hp2Tisa(const double& Hp_m) {
    if (Hp_m <= _Hp1_m) {
        return _T0_degK + _betaT * Hp_m; // standard temperature
    }
    else {
        return _T0_degK + env::atm::_betaT * env::atm::_Hp1_m; // standard temperature
    }
}
/* returns standard temperature [degK] based on pressure altitude */

double env::atm::Tisa2Hp(const double& Tisa_degK) {
    double Hp_m = (Tisa_degK - _T0_degK) / _betaT; // pressure altitude
    if (Hp_m < _Hp1_m) { // less than INSTEAD of less or equal than
        return Hp_m;
    }
    else {
        throw std::runtime_error("Pressure altitude can not be obtained from standard temperature above tropopause.");
    }
}
/* returns pressure altitude [m] based on standard temperature (no solution in stratosphere) */

double env::atm::Hp2T(const double& Hp_m, const double& DeltaT_degK) {
    return env::atm::Hp2Tisa(Hp_m) + DeltaT_degK; // temperature
}
/* returns temperature [degK] based on pressure altitude and temperature offset */

double env::atm::T2Hp(const double& T_degK, const double& DeltaT_degK) {
    double Hp_m = (T_degK - _T0_degK - DeltaT_degK) / _betaT; // pressure altitude
    if (Hp_m < _Hp1_m) { // less than INSTEAD of less or equal than
        return Hp_m;
    }
    else {
        throw std::runtime_error("Pressure altitude can not be obtained from temperature and temperature offset above tropopause.");
    }
}
/* returns pressure altitude [m] based on temperature and temperature offset (no solution in stratosphere) */

double env::atm::Hp2p(const double& Hp_m) {
    if (Hp_m <= _Hp1_m) {
        return _p0_pa * pow((1. + _betaT * Hp_m / _T0_degK), (_g0_R / _betaT)); // pressure
    }
    else {
        double p1_pa = env::atm::Hp2p(_Hp1_m); // pressure at tropopause
        return p1_pa * exp(_g0_R * (Hp_m - _Hp1_m) / env::atm::Hp2Tisa(_Hp1_m)); // pressure
    }
}
/* returns pressure [pa] based on pressure altitude. */

double env::atm::p2Hp(const double& p_pa) {
    double Hp_m = (_T0_degK / _betaT) * (pow(p_pa / _p0_pa, _betaT / _g0_R) - 1.); // pressure altitude below tropopause
    if (Hp_m <= _Hp1_m) {
        return Hp_m;
    }
    else {
        double p1_pa = env::atm::Hp2p(_Hp1_m); // pressure at tropopause
        return _Hp1_m + env::atm::Hp2Tisa(_Hp1_m) / _g0_R * log(p_pa / p1_pa); // pressure altitude above tropopause
    }
}
/* returns pressure altitude [m] based on pressure */

double env::atm::Hp2H(const double& Hp_m, const double& DeltaT_degK, const double& Deltap_pa) {
    double p_MSL_pa = _p0_pa + Deltap_pa; // mean sea level pressure
    double Hp_MSL_m = _T0_degK / _betaT * (pow((p_MSL_pa / _p0_pa), (_betaT / env::atm::_g0_R)) - 1); // mean sea level pressure altitude
    if (Hp_m <= _Hp1_m) {
        return Hp_m - Hp_MSL_m + DeltaT_degK / _betaT * log(env::atm::Hp2Tisa(Hp_m) / env::atm::Hp2Tisa(Hp_MSL_m)); // geopotential altitude
    }
    else {
        double H1_m = _Hp1_m - Hp_MSL_m + DeltaT_degK / _betaT * log(env::atm::Hp2Tisa(_Hp1_m) / env::atm::Hp2Tisa(Hp_MSL_m)); // tropopause geopotential altitude
        return H1_m + (Hp_m - _Hp1_m) * env::atm::Hp2T(_Hp1_m, DeltaT_degK) / env::atm::Hp2Tisa(_Hp1_m); // geopotential altitude
    }
}
/* returns geopotential altitude [m] based on pressure altitude, temperature offset, and pressure offset */

double env::atm::H2Hp(const double& H_m, const double& DeltaT_degK, const double& Deltap_pa) {
    Hp2H_aux Func;
    std::vector<double> par(2);
    par[0] = DeltaT_degK;
    par[1] = Deltap_pa;
    return Func.find_zero_secant(par, H_m, H_m-5, H_m+5, 1e-11, 100);
}
/* returns pressure altitude [m] based on geopotential altitude, temperature offset, and pressure offset (requires iteration) */

double env::atm::T2a(const double& T_degK) {
    return sqrt(T_degK * _kappa * _R);
}
/* returns speed of sound [mps] based on temperature */

double env::atm::pT2rho(const double& T_degK, const double& p_pa) {
    return p_pa / (T_degK * _R);
}
/* returns density [kgm3] based on temperature and pressure */

double env::atm::Hp2dTisadHp(const double& Hp_m) {
    return env::atm::Hp2dTdHp(Hp_m);
}
/* returns differential of standard temperature with pressure altitude [degK / m] based on pressure altitude */

double env::atm::Hp2dTdHp(const double& Hp_m) {
    return (Hp_m <= _Hp1_m)
        ? _betaT
        : 0.;
}
/* returns differential of temperature with pressure altitude [degK / m] based on pressure altitude */

double env::atm::Hp2dpdHp(const double& Hp_m) {
    return _g0_R * env::atm::Hp2p(Hp_m) / env::atm::Hp2Tisa(Hp_m);
}
/* returns differential of pressure with pressure altitude [pa / m] based on pressure altitude */

double env::atm::Hp2dHdHp(const double& Hp_m, const double& DeltaT_degK) {
    return env::atm::Hp2T(Hp_m, DeltaT_degK) / env::atm::Hp2Tisa(Hp_m);
}
/* returns differential of geopotential altitude with pressure altitude [-] based on pressure altitude and temperature offset */

/* ===== ===== ===== Non Static Functions ===== ===== ===== */
/* ======================================================== */

double env::atm::Hp2T(const double& Hp_m) const {
	return env::atm::Hp2T(Hp_m, _DeltaT_degK);
}
/* returns temperature [degK] based on pressure altitude */

double env::atm::T2Hp(const double& T_degK) const {
    return env::atm::T2Hp(T_degK, _DeltaT_degK);
}
/* returns pressure altitude [m] based on temperature (no solution in stratosphere) */

double env::atm::Hp2H(const double& Hp_m) const {
    return env::atm::Hp2H(Hp_m, _DeltaT_degK, _Deltap_pa);
}
/* returns geopotential altitude [m] based on pressure altitude */

double env::atm::H2Hp(const double& H_m) const {
    return env::atm::H2Hp(H_m, _DeltaT_degK, _Deltap_pa);
}
/* returns pressure altitude [m] based on geopotential altitude (requires iteration) */

double env::atm::Hp2dHdHp(const double& Hp_m) const {
    return env::atm::Hp2dHdHp(Hp_m, _DeltaT_degK);
}
/* returns differential of geopotential altitude with pressure altitude [-] based on pressure altitude */

/* ===== ===== ===== Private Functions ===== ===== ===== */
/* ===================================================== */

double env::atm::func::find_zero_secant(const std::vector<double>& par,
                                        const double& res,
                                        double x0,
                                        double x1,
                                        const double& tol,
                                        const unsigned short& iter) {
    double f0, f1, f2, x2;
    try {f0	= exec(x0, par);}
    catch (...) {
        throw std::runtime_error("Function evaluation error.");
    }
    try {f1	= exec(x1, par);}
    catch (...) {
        throw std::runtime_error("Function evaluation error.");
    }

    int c = 0;
    while (fabs(f1-res) > tol) {
        if (++c == iter) {
            throw std::runtime_error("Maximum number of iterations reached.");
        }
        if (f1 == f0) {throw std::runtime_error("Unable to find a zero.");}
        x2	= x0 - (f0 - res) * (x1 - x0) / (f1 - f0);
        try {f2	= exec(x2, par);}
        catch (...) {
            throw std::runtime_error("Function evaluation error.");
        }
        x0	= x1;
        f0	= f1;
        x1	= x2;
        f1	= f2;
    }
    return x1;
}

double env::atm::Hp2H_aux::exec(const double& Hp, const std::vector<double>& par) {
    // par[0] == DeltaT_degK
    // par[1] == Deltap_pa
    return env::atm::Hp2H(Hp, par[0], par[1]);
}

double env::atm::constructor_aux_eval(const double& Tisa_degK_msl, const double& Tisa_degK, const double& H_m) const {
    return((Tisa_degK - Tisa_degK_msl) + _DeltaT_degK * (log(Tisa_degK / Tisa_degK_msl)) - H_m * _betaT);
}

env::atm::constructor_aux::constructor_aux(const atm* const am_this) : _am_this(am_this) {}

double env::atm::constructor_aux::exec(const double& Tisa_degK_msl, const std::vector<double>& par) {
    // par[0] == Tisa
    // par[1] == H
    _Tisa_degK_msl = Tisa_degK_msl;
    _Tisa_degK     = par[0];
    _H_m           = par[1];
    return _am_this->constructor_aux_eval(_Tisa_degK_msl, _Tisa_degK, _H_m);
}
