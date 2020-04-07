#include "env/atm.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

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

void paper_relationships(const std::string&);
void paper_differentials(const std::string&);
void atmosphere_validation();
void atmosphere_identification();

int main(int argc, char **argv) {

    std::string st_output = "/data/PHD/outputs/"; // folder where text files are to be stored

    paper_relationships(st_output);
    paper_differentials(st_output);
    atmosphere_validation();
    atmosphere_identification();

	return 0;
}

void paper_relationships(const std::string& st_output) {
    std::string st_output_Tisa_Hp             = st_output + "atm/Tisa_Hp.txt";
    std::string st_output_T_Hp__DeltaT        = st_output + "atm/T_Hp__DeltaT.txt";
    std::string st_output_p_Hp                = st_output + "atm/p_Hp.txt";
    std::string st_output_H_Hp__DeltaT        = st_output + "atm/H_Hp__DeltaT.txt";
    std::string st_output_H_Hp__Deltap        = st_output + "atm/H_Hp__Deltap.txt";
    std::string st_output_H_Hp__DeltaT_Deltap = st_output + "atm/H_Hp__DeltaT_Deltap.txt";

    std::ofstream Oout_Tisa_Hp, Oout_T_Hp__DeltaT, Oout_p_Hp, Oout_H_Hp__DeltaT, Oout_H_Hp__Deltap, Oout_H_Hp__DeltaT_Deltap;
    Oout_Tisa_Hp.open(st_output_Tisa_Hp);
    Oout_T_Hp__DeltaT.open(st_output_T_Hp__DeltaT);
    Oout_p_Hp.open(st_output_p_Hp);
    Oout_H_Hp__DeltaT.open(st_output_H_Hp__DeltaT);
    Oout_H_Hp__Deltap.open(st_output_H_Hp__Deltap);
    Oout_H_Hp__DeltaT_Deltap.open(st_output_H_Hp__DeltaT_Deltap);

    env::atm Oam_n20(-20.0, 0.), Oam_n10(-10.0, 0.), Oam_p00(0., 0.), Oam_p10(10.0, 0.), Oam_p20(20.0, 0.);
    env::atm Oam_n5000(0,-5000.0), Oam_n2500(0,-2500.0), Oam_p2500(0,2500.0), Oam_p5000(0,5000.0);
    env::atm Oam_n20_n5000(-20.0,-5000.0), Oam_n20_p5000(-20.0,5000.0), Oam_p20_n5000(20.0,-5000.0), Oam_p20_p5000(20.0,5000.0);

    int nel = 151;
    std::vector<double>VHp_m(nel), VTisa_degK(nel), Vp_pa(nel);
    std::vector<double>VT_degK_n20(nel), VT_degK_n10(nel), VT_degK_p00(nel), VT_degK_p10(nel), VT_degK_p20(nel);
    std::vector<double>VH_m_n20(nel), VH_m_n10(nel), VH_m_p00(nel), VH_m_p10(nel), VH_m_p20(nel);
    std::vector<double>VH_m_n5000(nel), VH_m_n2500(nel), VH_m_p2500(nel), VH_m_p5000(nel);
    std::vector<double>VH_m_n20_n5000(nel), VH_m_n20_p5000(nel), VH_m_p20_n5000(nel), VH_m_p20_p5000(nel);

    for (int i = 0; i != nel; ++i) {
        VHp_m[i] = 100.0 * i;
        VTisa_degK[i]  = env::atm::Hp2Tisa(VHp_m[i]);
        VT_degK_n20[i] = Oam_n20.Hp2T(VHp_m[i]);
        VT_degK_n10[i] = Oam_n10.Hp2T(VHp_m[i]);
        VT_degK_p00[i] = Oam_p00.Hp2T(VHp_m[i]);
        VT_degK_p10[i] = Oam_p10.Hp2T(VHp_m[i]);
        VT_degK_p20[i] = Oam_p20.Hp2T(VHp_m[i]);
        Vp_pa[i]       = env::atm::Hp2p(VHp_m[i]);
        VH_m_n20[i]    = Oam_n20.Hp2H(VHp_m[i]);
        VH_m_n10[i]    = Oam_n10.Hp2H(VHp_m[i]);
        VH_m_p00[i]    = Oam_p00.Hp2H(VHp_m[i]);
        VH_m_p10[i]    = Oam_p10.Hp2H(VHp_m[i]);
        VH_m_p20[i]    = Oam_p20.Hp2H(VHp_m[i]);
        VH_m_n5000[i]  = Oam_n5000.Hp2H(VHp_m[i]);
        VH_m_n2500[i]  = Oam_n2500.Hp2H(VHp_m[i]);
        VH_m_p2500[i]  = Oam_p2500.Hp2H(VHp_m[i]);
        VH_m_p5000[i]  = Oam_p5000.Hp2H(VHp_m[i]);
        VH_m_n20_n5000[i] = Oam_n20_n5000.Hp2H(VHp_m[i]);
        VH_m_n20_p5000[i] = Oam_n20_p5000.Hp2H(VHp_m[i]);
        VH_m_p20_n5000[i] = Oam_p20_n5000.Hp2H(VHp_m[i]);
        VH_m_p20_p5000[i] = Oam_p20_p5000.Hp2H(VHp_m[i]);

        Oout_Tisa_Hp << std::fixed << std::showpos
                     << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                     << std::setprecision(3) << std::setw(15) << VTisa_degK[i]
                     << std::endl;

        Oout_T_Hp__DeltaT << std::fixed << std::showpos
                          << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                          << std::setprecision(3) << std::setw(15) << VT_degK_n20[i]
                          << std::setprecision(3) << std::setw(15) << VT_degK_n10[i]
                          << std::setprecision(3) << std::setw(15) << VT_degK_p00[i]
                          << std::setprecision(3) << std::setw(15) << VT_degK_p10[i]
                          << std::setprecision(3) << std::setw(15) << VT_degK_p20[i]
                          << std::endl;

        Oout_p_Hp << std::fixed << std::showpos
                  << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                  << std::setprecision(5) << std::setw(15) << Vp_pa[i] * 1e-3
                  << std::endl;

        Oout_H_Hp__DeltaT << std::fixed << std::showpos
                          << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_n20[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_n10[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_p00[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_p10[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_p20[i] * 1e-3
                          << std::endl;

        Oout_H_Hp__Deltap << std::fixed << std::showpos
                          << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_n5000[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_n2500[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_p00[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_p2500[i] * 1e-3
                          << std::setprecision(3) << std::setw(12) << VH_m_p5000[i] * 1e-3
                          << std::endl;

        Oout_H_Hp__DeltaT_Deltap << std::fixed << std::showpos
                                 << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                                 << std::setprecision(3) << std::setw(12) << VH_m_n20_n5000[i] * 1e-3
                                 << std::setprecision(3) << std::setw(12) << VH_m_n20_p5000[i] * 1e-3
                                 << std::setprecision(3) << std::setw(12) << VH_m_p00[i] * 1e-3
                                 << std::setprecision(3) << std::setw(12) << VH_m_p20_n5000[i] * 1e-3
                                 << std::setprecision(3) << std::setw(12) << VH_m_p20_p5000[i] * 1e-3
                                 << std::endl;
    }
    Oout_Tisa_Hp.close();
    Oout_T_Hp__DeltaT.close();
    Oout_p_Hp.close();
    Oout_H_Hp__DeltaT.close();
    Oout_H_Hp__Deltap.close();
    Oout_H_Hp__DeltaT_Deltap.close();

    std::cout << st_output_Tisa_Hp << std::endl;

}

void paper_differentials(const std::string& st_output) {
    std::string st_output_dTisadHp_Hp      = st_output + "atm/dTisadHp_Hp.txt";
    std::string st_output_dTdHp_Hp         = st_output + "atm/dTdHp_Hp.txt";
    std::string st_output_dpdHp_Hp         = st_output + "atm/dpdHp_Hp.txt";
    std::string st_output_dHdHp_Hp__DeltaT = st_output + "atm/dHdHp_Hp__DeltaT.txt";

    std::ofstream Oout_dTisadHp_Hp, Oout_dTdHp_Hp, Oout_dpdHp_Hp, Oout_dHdHp_Hp__DeltaT;
    Oout_dTisadHp_Hp.open(st_output_dTisadHp_Hp);
    Oout_dTdHp_Hp.open(st_output_dTdHp_Hp);
    Oout_dpdHp_Hp.open(st_output_dpdHp_Hp);
    Oout_dHdHp_Hp__DeltaT.open(st_output_dHdHp_Hp__DeltaT);

    env::atm Oam_n20(-20.0, 0.), Oam_n10(-10.0, 0.), Oam_p00(0., 0.), Oam_p10(10.0, 0.), Oam_p20(20.0, 0.);

    int nel = 151;
    std::vector<double>VHp_m(nel), VdTisadHp_degKm(nel), VdTdHp_degKm(nel), VdpdHp_pam(nel);
    std::vector<double>VdHdHp_n20(nel), VdHdHp_n10(nel), VdHdHp_p00(nel), VdHdHp_p10(nel), VdHdHp_p20(nel);

    for (int i = 0; i != nel; ++i) {
        VHp_m[i] = 100.0 * i;
        VdTisadHp_degKm[i] = env::atm::Hp2dTisadHp(VHp_m[i]);
        VdTdHp_degKm[i]    = env::atm::Hp2dTdHp(VHp_m[i]);
        VdpdHp_pam[i]      = env::atm::Hp2dpdHp(VHp_m[i]);
        VdHdHp_n20[i]      = Oam_n20.Hp2dHdHp(VHp_m[i]);
        VdHdHp_n10[i]      = Oam_n10.Hp2dHdHp(VHp_m[i]);
        VdHdHp_p00[i]      = Oam_p00.Hp2dHdHp(VHp_m[i]);
        VdHdHp_p10[i]      = Oam_p10.Hp2dHdHp(VHp_m[i]);
        VdHdHp_p20[i]      = Oam_p20.Hp2dHdHp(VHp_m[i]);

        Oout_dTisadHp_Hp << std::fixed << std::showpos
                         << std::setprecision(3) << std::setw(10) << VHp_m[i] * 1e-3
                         << std::setprecision(1) << std::setw(15) << VdTisadHp_degKm[i] * 1e3
                         << std::endl;

        Oout_dTdHp_Hp << std::fixed << std::showpos
                      << std::setprecision(3) << std::setw(10) << VHp_m[i] * 1e-3
                      << std::setprecision(1) << std::setw(15) << VdTdHp_degKm[i] * 1e3
                      << std::endl;

        Oout_dpdHp_Hp << std::fixed << std::showpos
                      << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                      << std::setprecision(3) << std::setw(15) << VdpdHp_pam[i]
                      << std::endl;

        Oout_dHdHp_Hp__DeltaT << std::fixed << std::showpos
                              << std::setprecision(1) << std::setw(10) << VHp_m[i] * 1e-3
                              << std::setprecision(5) << std::setw(15) << VdHdHp_n20[i]
                              << std::setprecision(5) << std::setw(15) << VdHdHp_n10[i]
                              << std::setprecision(5) << std::setw(15) << VdHdHp_p00[i]
                              << std::setprecision(5) << std::setw(15) << VdHdHp_p10[i]
                              << std::setprecision(5) << std::setw(15) << VdHdHp_p20[i]
                              << std::endl;

        if (fabs(VHp_m[i] - 11000) < 1e-4) {
            Oout_dTisadHp_Hp << std::fixed << std::showpos
                             << std::setprecision(3) << std::setw(10) << VHp_m[i] * 1e-3 + 1e-3
                             << std::setprecision(1) << std::setw(15) << VdTisadHp_degKm[i+1] * 1e3
                             << std::endl;

            Oout_dTdHp_Hp << std::fixed << std::showpos
                          << std::setprecision(3) << std::setw(10) << VHp_m[i] * 1e-3 + 1e-3
                          << std::setprecision(1) << std::setw(15) << VdTdHp_degKm[i+1] * 1e3
                          << std::endl;
        }


    }
    Oout_dTisadHp_Hp.close();
    Oout_dTdHp_Hp.close();
    Oout_dpdHp_Hp.close();
    Oout_dHdHp_Hp__DeltaT.close();

    std::cout << st_output_dTdHp_Hp << std::endl;

}

void atmosphere_validation() {
    double DeltaT_degK = 12.5;
    double Deltap_pa   = 3800.0;
    env::atm Oam(DeltaT_degK, Deltap_pa);

    double Hp_m_trop = 5000.0;

    double Tisa_trop_degK = Oam.Hp2Tisa(Hp_m_trop);
    double T_trop_degK    = Oam.Hp2T(Hp_m_trop);
    double p_trop_pa      = Oam.Hp2p(Hp_m_trop);
    double H_trop_m       = Oam.Hp2H(Hp_m_trop);
    double rho_trop_kgm3  = Oam.pT2rho(T_trop_degK, p_trop_pa);
    double a_mps          = Oam.T2a(T_trop_degK);

    double Hp_m_trop_2    = Oam.Tisa2Hp(Tisa_trop_degK);
    double Hp_m_trop_3    = Oam.T2Hp(T_trop_degK);
    double Hp_m_trop_4    = Oam.p2Hp(p_trop_pa);
    double Hp_m_trop_5    = Oam.H2Hp(H_trop_m);

    std::cout << Hp_m_trop - Hp_m_trop_2 << std::endl;
    std::cout << Hp_m_trop - Hp_m_trop_3 << std::endl;
    std::cout << Hp_m_trop - Hp_m_trop_4 << std::endl;
    std::cout << Hp_m_trop - Hp_m_trop_5 << std::endl;

    double Hp_m_strat = 13500.0;

    double Tisa_strat_degK = Oam.Hp2Tisa(Hp_m_strat);
    double T_strat_degK    = Oam.Hp2T(Hp_m_strat);
    double p_strat_pa      = Oam.Hp2p(Hp_m_strat);
    double H_strat_m       = Oam.Hp2H(Hp_m_strat);

    double Hp_m_strat_4    = Oam.p2Hp(p_strat_pa);
    double Hp_m_strat_5    = Oam.H2Hp(H_strat_m);

    std::cout << Hp_m_strat - Hp_m_strat_4 << std::endl;
    std::cout << Hp_m_strat - Hp_m_strat_5 << std::endl;
}

void atmosphere_identification() {
    double DeltaT_degK = 12.5;
    double Deltap_pa   = 3800.0;
    env::atm Oam(DeltaT_degK, Deltap_pa);

    double Hp_m = 1500.0;
    double T_degK = Oam.Hp2T(Hp_m);
    double p_pa   = Oam.Hp2p(Hp_m);
    double H_m    = Oam.Hp2H(Hp_m);

    // identify new atmosphere based on these parameters
    env::atm Oam2(T_degK, p_pa, H_m);

    double DeltaT_degK_2 = Oam2.get_DeltaT_degK();
    double Deltap_pa_2   = Oam2.get_Deltap_pa();

    std::cout << DeltaT_degK - DeltaT_degK_2 << std::endl;
    std::cout << Deltap_pa - Deltap_pa_2 << std::endl;
}