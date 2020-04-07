#ifndef ENV_ATM
#define ENV_ATM

#include "env.h"
#include <vector>

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

namespace env {

// CLASS ATM
// =========
// =========

class ENV_API atm {
public:
    /**< ===== ===== ===== Atmospheric constants ===== ===== ===== */
	/**< air specific constant [J/(kg*degK)] = [m2/(degK*s2)] */
	static const double _R;
	/**< air adiabatic index [-] */
	static const double _kappa;
	/**< troposphere thermal gradient [degK/m] */
	static const double _betaT;
    /**< standard acceleration of free fall [m/s2] */
    static const double _g0_mps2;
    /**< division of g0 by R, with negative sign [degK/m] */
    static const double	_g0_R;
    /**< tropopause pressure altitude */
    static const double _Hp1_m;

    /**< ===== ===== ===== Standard mean sea level conditions ===== ===== ==== */
    /**< Conditions at standard mean sea level Hp = 0 if both temperature and pressure offsets are zero */
	/**< standard mean sea level pressure */
	static const double _p0_pa;
	/**< standard mean sea level temperature */
	static const double _T0_degK;
private:
	/**< temperature offset */
	double _DeltaT_degK;
	/**< pressure offset */
	double _Deltap_pa;
public:
    /**< ===== ===== ===== Constructors ===== ===== ===== */
	/**< constructor based on temperature and pressure differentials at mean sea level. */
	explicit atm(const double& DeltaT_degK = 0, const double& Deltap_pa = 0);
    /**< constructor based on given point temperature, pressure, and geopotential altitude. This point, normally at
     * airport or meteorological station, is assumed to be in the troposphere). */
	atm(const double& T_degK_datum, const double& p_pa_datum, const double& H_m_datum);

    /**< copy constructor */
	atm(const atm&) = default;
    /**< move constructor */
    atm(atm&&) = default;
    /**< destructor */
    ~atm() = default;
    /**< copy assignment */
    atm& operator=(const atm&) = default;
    /**< move assignment */
    atm& operator=(atm&&) = default;

    /**< ===== ===== ===== Get & Set ===== ===== ===== */
	/**< sets temperature and pressure offsets at mean sea level */
    void set(const double& DeltaT_degK, const double& Deltap_pa);
	/**< returns value of the temperature offset at mean sea level */
    const double& get_DeltaT_degK() const {return _DeltaT_degK;}
	/**< returns value of the pressure offset at mean sea level. */
    const double& get_Deltap_pa() const {return _Deltap_pa;}

    /**< ===== ===== ===== Static Functions ===== ===== ===== */
    /**< ==================================================== */
    /**< returns standard temperature [degK] based on pressure altitude */
    static double Hp2Tisa(const double& Hp_m);
    /**< returns pressure altitude [m] based on standard temperature (no solution in stratosphere) */
    static double Tisa2Hp(const double& Tisa_degK);

    /**< returns temperature [degK] based on pressure altitude and temperature offset */
    static double Hp2T(const double& Hp_m, const double& DeltaT_degK);
    /**< returns pressure altitude [m] based on temperature and temperature offset (no solution in stratosphere) */
    static double T2Hp(const double& T_degK, const double & DeltaT_degK);

	/**< returns pressure [pa] based on pressure altitude. */
    static double Hp2p(const double& Hp_m);
	/**< returns pressure altitude [m] based on pressure */
    static double p2Hp(const double& p_pa);

    /**< returns geopotential altitude [m] based on pressure altitude, temperature offset, and pressure offset */
    static double Hp2H(const double& Hp_m, const double& DeltaT_degK, const double& Deltap_pa);
    /**< returns pressure altitude [m] based on geopotential altitude, temperature offset, and pressure offset (requires iteration) */
    static double H2Hp(const double& H_m, const double& DeltaT_degK, const double& Deltap_pa);

    /**< returns speed of sound [mps] based on temperature */
    static double T2a(const double& T_degK);
    /**< returns density [kgm3] based on temperature and pressure */
    static double pT2rho(const double& T_degK, const double& p_pa);

    /**< returns differential of standard temperature with pressure altitude [degK / m] based on pressure altitude */
    static double Hp2dTisadHp(const double& Hp_m);
    /**< returns differential of temperature with pressure altitude [degK / m] based on pressure altitude */
    static double Hp2dTdHp(const double& Hp_m);
    /**< returns differential of pressure with pressure altitude [pa / m] based on pressure altitude */
    static double Hp2dpdHp(const double& Hp_m);
    /**< returns differential of geopotential altitude with pressure altitude [-] based on pressure altitude and temperature offset */
    static double Hp2dHdHp(const double& Hp_m, const double& DeltaT_degK);

    /**< ===== ===== ===== Non Static Functions ===== ===== ===== */
    /**< ======================================================== */
    /**< returns temperature [degK] based on pressure altitude */
    double Hp2T(const double& Hp_m) const;
    /**< returns pressure altitude [m] based on temperature (no solution in stratosphere) */
    double T2Hp(const double& T_degK) const;

    /**< returns geopotential altitude [m] based on the pressure altitude */
    double Hp2H(const double& Hp_m) const;
    /**< returns pressure altitude [m] based on geopotential altitude (requires iteration) */
    double H2Hp(const double& H_m) const;

    /**< returns differential of geopotential altitude with pressure altitude [-] based on pressure altitude */
    double Hp2dHdHp(const double& Hp_m) const;
private:

    // class intended to obtain the zero of a given function by means of the secant method */
    class func {
    public:
        virtual double exec	(const double&, const std::vector<double>&) = 0;
        virtual double find_zero_secant(const std::vector<double>& par, const double& res, double x1, double x2, const double& tol, const unsigned short& iter);
    };

    /**< class intended to assist in iteration */
    class Hp2H_aux : public func {
    public:
        double exec(const double& Hp, const std::vector<double>& par) override;
    };

    /**< class intended to assist in iteration */
    double constructor_aux_eval(const double& Tisa_degK_msl, const double& Tisa_degK, const double& H_m) const;
    class constructor_aux : public func {
    private:
        const atm* const _am_this;
        double _Tisa_degK_msl;
        double _Tisa_degK;
        double _H_m;
    public:
        explicit constructor_aux(const atm* const);
        double exec(const double& Tisa_degK_msl, const std::vector<double>& par) override;
    };
}; // closes class atm

} // closes namespace env


#endif


