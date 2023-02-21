#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

//#define NUM_THREADS 4

#define dt   0.0025   //0.05   // time step (ms)
#define dx   0.04     // mm

#define NUM_CELLS   147609
double h_sparse[NUM_CELLS]; // h, Fast inactivation gating variable for I_Na (dimensionless)
double m_sparse[NUM_CELLS]; // m, Activation gating variable for I_Na (dimensionless)
double j_sparse[NUM_CELLS]; // j, Slow inactivation gating variable for I_Na (dimensionless)
double oa_sparse[NUM_CELLS]; // Activation gating variable for I_to (dimensionless)
double oi_sparse[NUM_CELLS]; // Inactivation gating variable for I_to (dimensionless)
double ua_sparse[NUM_CELLS]; // Activation gating variable for I_Kur (dimensionless)
double ui_sparse[NUM_CELLS]; // Inactivation gating variable for I_Kur (dimensionless)
double xr_sparse[NUM_CELLS]; // xr, Rapidly activating K time-dependant activation gate (dimensionless)
double paF_sparse[NUM_CELLS];
double paS_sparse[NUM_CELLS];
double xs_sparse[NUM_CELLS]; // xs, Activation gating variable for I_Ks (dimensionless)
double d_sparse[NUM_CELLS]; // d, Activation gating variable for I_CaL  (dimensionless)
double f_sparse[NUM_CELLS]; // f, Voltage-dependent inactivation gating variable for I_CaL (dimensionless)
double fCa_sparse[NUM_CELLS]; // fCa, Ca-dependent inactivation gating variable for I_CaL (dimensionless)
double dT_sparse[NUM_CELLS];
double fT_sparse[NUM_CELLS];
double y_sparse[NUM_CELLS];
double a_sparse[NUM_CELLS];
double Na_i_sparse[NUM_CELLS]; // Intracellular concentration of Na (mM)
double K_i_sparse[NUM_CELLS]; // Intracellular concentration of K (mM)

double RyR_a_NJ_sparse[NUM_CELLS];
double RyR_o_NJ_sparse[NUM_CELLS];
double RyR_c_NJ_sparse[NUM_CELLS];
double Ca_i_NJ_sparse[NUM_CELLS];
double Ca_sr_NJ_sparse[NUM_CELLS];
double Ca_serca_NJ_sparse[NUM_CELLS];
double RyR_a_SS_sparse[NUM_CELLS];
double RyR_o_SS_sparse[NUM_CELLS];
double RyR_c_SS_sparse[NUM_CELLS];
double Ca_i_SS_sparse[NUM_CELLS];
double Ca_sr_SS_sparse[NUM_CELLS];
double Ca_serca_SS_sparse[NUM_CELLS];
double fCQ_sparse[NUM_CELLS];
double fTC_sparse[NUM_CELLS];


double D1m_sparse[NUM_CELLS];
double D2m_sparse[NUM_CELLS];
double V_sparse[NUM_CELLS]; // Transmembrane potential (mV), timestep n
double V_new_sparse[NUM_CELLS]; // Transmembrane potential (mV), timestep n+1
double fx_sparse[NUM_CELLS];
double fy_sparse[NUM_CELLS];
double fz_sparse[NUM_CELLS];

short int geo_x_index[NUM_CELLS];
short int geo_y_index[NUM_CELLS];
short int geo_z_index[NUM_CELLS];
short int cell_map_sparse[NUM_CELLS];
char IC_map[NUM_CELLS];
int geo_adj[NUM_CELLS][18];

double D1_cube[3][3][3]; 
double D2_cube[3][3][3];
double fx_cube[3][3][3];
double fy_cube[3][3][3];
double fz_cube[3][3][3];
double  V_cube[3][3][3];
double   I_Na; // Fast inward Na current (pA)
double   I_K1; // Inward rectifier K current (pA)
double   I_to; // Transient outward K current (pA)
double  I_Kur; // Ultrarapid delayed rectifier K current (pA)
double   I_Kr; // Rapid delayed rectifier K current (pA)
double   I_Ks; // Slow delayed rectifier K current (pA)
double  I_CaL; // L-type inward Ca current (pA)
double  I_CaT;
double  I_f;
double  I_NaK; // Na-K pump current (pA)
double I_NaCa; // Na/Ca exchanger current (pA)
double  I_bCa; // Background Ca current (pA)
double  I_bNa; // Background Na current (pA)
double  I_pCa; // Sarcoplasmic Ca pump current (pA)
double I_KACh; // Acetylcholine-activated K current (pA)
double  I_ion; // Total ionic current (pA)

double D1_atria_var;
double D2_atria_var;
double D1_SAN_var;
double D2_SAN_var;
double D1_var;
double D2_var;

double const R = 8.3143; // Gas constant (J/mol*K)
double const F = 96.4867; // Faraday constant (C/mmol)
double const T = 310.0; // Temperature (K)

/* Cell Geometry */
double C_m;         // Membrane capacitance (pF)
double const Volume_i = 13668.0;    // Intracellular volume (um^3)
double const L_cell = 67.0;   // micrometre (in Cell_parameters)
double const L_sub = 0.02;   // micrometre (in Cell_parameters)
double const R_cell = 3.9;   // micrometre (in Cell_parameters)
double const V_i_part = 0.46;   // dimensionless (in Cell_parameters)
double const V_jsr_part = 0.0012;   // dimensionless (in Cell_parameters)
double const V_nsr_part = 0.0116;   // dimensionless (in Cell_parameters)
double V_sub;
double V_cell;
double V_nsr;
double V_i;
double V_jsr;

/* Ion Valences */
double const z_Na = 1.0; // Na valence
double const z_K = 1.0;  // K valence
double const z_Ca = 2.0; // Ca valence

/* Ion Concentrations */
double const K_o = 5.4;      // Extracellular concentration of K (mM)
double const K_i_SAN = 140.0;
double const Na_o = 140.0;     // Extracellular concentration of Na (mM)
double const Na_i_SAN = 5.0;
double const Ca_o = 1.8;     // Extracellular concentration of Ca (mM)

/* Ion Equilibrium Potentials */
double E_K;      // Equilibrium potential for K (mV)
double E_Na;     // Equilibrium potential for Na (mV)
double E_Ca;     // Equilibrium potential for Ca (mV)
double E_Na_SAN;
double E_mh_SAN;     
double E_K_SAN;     
double E_Ks_SAN;     

double const SLlow = 165;
double const KdSLlow = 1.1;
double const SLhigh = 13;
double const KdSLhigh = 13e-3;
double const BCa = 24e-3;
double const KdBCa = 2.38e-3;
double const CSQN = 6.7;
double const KdCSQN = 0.8;
double const DCa = 780;
double const ASS_NJ = 2492.324412;
double const XSS_NJ = 0.822500;
double const k3 = 2.314815;    
double const k4 = 7.5;            // % pump rate
double const k1 = 7500000.000000; 
double const k2 = 0.468750;       
double const Cpumps = 40e-3;
double const V_NJ = 6 * 0.002531;
double const V_SS = 2 * 4.99232e-5;
double const kSRleak = 6e-3;
double const TRyRadapt = 1 * 0.25;
double const TRyRactSS = 0.005;
double const TRyRactNJ = 5e-3;
double const TRyRinactSS = 15e-3;
double const TRyRinactNJ = 2 * 15e-3;
double const DCaSR = 44;
double const r = 1.625;
double const V_SR_NJ = 2 * 0.000057;
double const V_SR_SS = 2 * 0.000080;

double BSS;
double B_i_NJ;
double B_sr_NJ;
double B_sr_SS;
double I_SS_NJ;
double I_serca_SR_NJ;
double I_serca_bulk_NJ;
double I_serca_SR_SS;
double I_serca_bulk_SS;
double RyR_serca_SS;
double RyR_a_SS_infinity;
double RyR_o_SS_infinity;
double RyR_c_SS_infinity;
double I_rel_SS;
double RyR_serca_NJ;
double RyR_a_NJ_infinity;
double RyR_o_NJ_infinity;
double RyR_c_NJ_infinity;
double I_rel_NJ;
double I_SRleak_NJ;
double I_SRleak_SS;
double I_Ca_NJ;
double I_Ca_SS;
double I_SRCa_NJ;
double I_SRCa_SS;
double dCa_SR_NJ;
double dCa_SR_SS;

/* Ca system SAN */
double j_SRCarel, diff, kCaSR, koSRCa, kiSRCa, delta_R, delta_O, delta_I, delta_RI, P_tot;
double delta_fTC, delta_fTMC, delta_fTMM , delta_fCMi , delta_fCMs , delta_fCQ;
double j_Ca_dif, ACh_block_P_up, j_up, delta_Ca_i, delta_Ca_sub, j_tr, delta_Ca_nsr, delta_Ca_jsr;


/* Fast inward Sodium/Na current */
double g_Na;       // Maximal I_Na conductance (nS/pF)
double alpha_h;    // Forward rate constant for gating variable h in I_Na (1/ms)
double beta_h;     // Backward rate constant for gating variable h in I_Na (1/ms)
double alpha_j;    // Forward rate constant for gating variable j in I_Na (1/ms)
double beta_j;     // Backward rate constant for gating variable j in I_Na (1/ms)
double alpha_m;    // Forward rate constant for gating variable m in I_Na (1/ms)
double beta_m;     // Backward rate constant for gating variable m in I_Na (1/ms)
double m_infinity; // Steady-state relation for gating variable m in I_Na (dimensionless)
double h_infinity; // Steady-state relation for gating variable h in I_Na (dimensionless)
double j_infinity; // Steady-state relation for gating variable j in I_Na (dimensionless)
double tau_m;      // Time constant for gating variable m in I_Na (ms)
double tau_h;      // Time constant for gating variable h in I_Na (ms)
double tau_j;      // Time constant for gating variable j in I_Na (ms)
double g_Na_SAN;

/* Inward rectifier K current */
double g_K1; // Maximal I_K1 conductance (nS/pF)
double g_K1_SAN;

/* Transient Outward K Current */
double g_to;        // Maximal I_to conductance (nS/pF)
double g_to_SAN; 
double tau_oa;      // Time constant for gating variable oa in I_to (ms)
double oa_infinity; // Steady-state relation for gating variable oa in I_to (dimensionless)

double tau_oi;      // Time constant for gating variable oi in I_to (ms)
double oi_infinity; // Steady-state relation for gating variable oi in I_to (dimensionless)
double const K_Q10 = 3.0;       // Temprature scaling factor for I_Kur and I_to kinetics (dimensionless)

/* Ultrarapid delayed rectifier K current */
double g_Kur;       // Maximal I_Kur conductance (nS/pF)
double g_Kur_SAN;
double ua_infinity; // Steady-state relation for gating variable ua in I_Kur (dimensionless)
double tau_ua;      // Time constant for gating variable ua in I_Kur (ms)
double ui_infinity; // Steady-state value of activation gate ui
double tau_ui;      // Steady-state relation for gating variable ui in I_Kur (dimensionless)

/* Rapid delayed rectifier K current */
double g_Kr;        // Maximal I_Kr conductance (nS/pF)
double xr_infinity; // Steady-state relation for gating variable xr in I_Kr (dimensionless)
double tau_xr;      // Time constant for gating variable xr in I_Kr (ms)
double alpha_xr;    // Forward rate constant for gating variable xr in I_Kr (1/ms)
double beta_xr;     // Backward rate constant for gating variable xr in I_Kr (1/ms)
double g_Kr_SAN;
double pa_infinity;
double tau_paF;
double tau_paS;

/* Slow delayed rectifier K current */
double g_Ks;        // Maximal I_Ks conductance (nS/pF)
double g_Ks_SAN;
double xs_infinity; // Steady-state relation for gating variable xs in I_Ks (dimensionless)
double tau_xs;      // Time constant for gating variable xs in I_Ks (ms)
double alpha_xs;    // Forward rate constant for gating variable xs in I_Ks (1/ms)
double beta_xs;     // Backward rate constant for gating variable xs in I_Ks (1/ms)

/* L-type inward Ca current */
double alpha_d, beta_d;
double tau_d;        // Time constant for gating variable d in I_CaL (ms)
double tau_f;        // Time constant for gating variable f in I_CaL (ms)
double tau_fCa;      // Time constant for gating variable fCa in I_CaL (ms)
double d_infinity;   // Steady-state relation for gating variable d in I_CaL (dimensionless)
double f_infinity;   // Steady-state relation for gating variable f in I_CaL (dimensionless)
double fCa_infinity; // Steady-state relation for gating variable fCa in I_CaL (dimensionless)
double g_CaL;        // Maximal I_CaL conductance (nS/pF)
double P_CaL;
double I_siCa, I_siNa, I_siK;
double adVm, bdVm;
double ACh_block_I_CaL;

/* I_CaT */
double dT_infinity, fT_infinity;
double tau_dT, tau_fT;
double P_CaT;

/* funny current */
double const Km_f = 45.0;   // millimolar (in i_f)
double const alpha = 0.5927;   // dimensionless (in i_f)
double const g_f = 0.00427;   // microS (in i_f)
double G_f, G_f_K, G_f_Na, g_f_Na, g_f_K, I_f_Na, I_f_K;
double tau_y, y_infinity;
double ACh_shift_I_f;

/* Na-K pump current */
double f_NaK;    // Voltage-dependance parameter of I_NaK (dimensionless)
double sigma;    // Na_o dependance parameter for f_NaK (dimensionless)
double I_NaKmax; // Maximal I_NaK (pA/pF)
double I_NaKmax_SAN;
double const K_mNai = 0.95 * 10.0;   // Na_i half-saturation constant for I_NaK (mM)
double const K_mKo = 1.5;    // K_o half-saturation constant for I_NaK (mM)
double const K_mKo_SAN = 1.4;
double const K_mNai_SAN = 14.0;

/* Na/Ca exchanger current */
double const K_mNa = 87.5;     // Na_o saturation constant for I_NaCa (mM)
double const k_sat = 0.1;     // Saturation factor for I_NaCa (dimensionless)
double const K_mCa = 1.38;     // Ca saturation factor for NaCa exchanger (mM)
double const Gamma = 0.35;     // Voltage-dependence parameter for I_NaCa (dimensionless)
double I_NaCamax; // Maximal I_NaCa (pA/pF)
double I_NaCamax_SAN;
double di, dodo;
double k12, k14, k21, k23, k32, k34, k41, k43;
double x1, x2, x3, x4;


/* Background Ca current */
double g_bCa; // Maximal I_bCa conductance (nS/pF)

/* Background Na current */
double g_bNa; // Maximal I_bNa conductance (nS/pF)

/* Sarcoplasmic Ca pump current */
double I_pCamax; // Maximal I_pCa (pA/pF)

/* Acetylcholine-Activated Potassium Current */
double ACh;        // Acetylcholine concentration (mM)
double alpha_a;
double beta_a;
double a_infinity;
double tau_a;
double g_KACh_SAN;

double D1, D2;
double fx, fy, fz;
double dfxdx, dfxdy, dfxdz, dfydx, dfydy, dfydz, dfzdx, dfzdy, dfzdz;
double D000, Dp00, Dm00, D0p0, D0m0, D00p, D00m, Dpp0, Dpm0, Dmp0, Dmm0, Dp0p, Dp0m, Dm0p, Dm0m, D0pp, D0pm, D0mp, D0mm;

void Courtemanche_Calculations (int sparse_index);
/* Ion Concentration Functions */
void Calculate_Na_i_Concentration       (int sparse_index); // Calculates new intracellular concentration of Na (mM)
void Calculate_K_i_Concentration        (int sparse_index); // Calculates new intracellular concentration of K (mM)
void Calculate_Ca_System_Concentrations (int sparse_index);
/* Ion Current Functions */     
void Calculate_I_Na   (int sparse_index); // Calculates fast inward Na current (pA)
void Calculate_I_K1   (int sparse_index); // Calculates inward rectifier K current (pA)
void Calculate_I_to   (int sparse_index); // Calculates transient outward K current (pA)
void Calculate_I_Kur  (int sparse_index); // Calculates ultrarapid delayed rectifier K current (pA)
void Calculate_I_Kr   (int sparse_index); // Calculates rapid delayed rectifier K current (pA)
void Calculate_I_Ks   (int sparse_index); // Calculates slow delayed rectifier K current (pA)
void Calculate_I_CaL  (int sparse_index); // Calculates L-type inward Ca current (pA)
void Calculate_I_CaT  (int sparse_index);
void Calculate_I_f    (int sparse_index);
void Calculate_I_NaK  (int sparse_index); // Calculates Na-K pump current (pA)
void Calculate_I_NaCa (int sparse_index); // Calculates Na/Ca exchanger current (pA)
void Calculate_I_bCa  (int sparse_index); // Calculates background Ca current (pA)
void Calculate_I_bNa  (int sparse_index); // Calculates background Na current (pA)
void Calculate_I_pCa  (int sparse_index); // Calculates sarcoplasmic Ca pump current (pA)
void Calculate_I_KACh (int sparse_index); // Calculates acetylcholine-activated K current (pA)
void Calculate_I_ion  (int sparse_index); // Calculates total ionic Current

#pragma omp threadprivate(D1, D2, D1_atria_var, D2_atria_var, D1_SAN_var, D2_SAN_var, D1_var, D2_var) 

void min_of_2(double D1_value_1, double D2_value_1, double D1_value_2, double D2_value_2){
	if(D1_value_1 < D1_value_2){
		D1 = D1_value_1;
		D2 = D2_value_1;
	}
	else{
		D1 = D1_value_2;
		D2 = D2_value_2;
	}
}

void min_of_3(double D1_value_1, double D2_value_1, double D1_value_2, double D2_value_2, double D1_value_3, double D2_value_3) {
	if((D1_value_1 < D1_value_2) && (D1_value_1 < D1_value_3)){
		D1 = D1_value_1;
		D2 = D2_value_1;
	}
	else{ 
		if(D1_value_2 < D1_value_3){
			D1 = D1_value_2; 
			D2 = D2_value_2;
		} 
		else{ 
			D1 = D1_value_3; 
			D2 = D2_value_3;
		}
	}
}

int disable_burst_pacing = 0; //default = 0
double pow_1p5(double base) { return base*sqrt(base); }
double pow_2(double base) { return base*base; }
double pow_3(double base) { return base*base*base; }
double norm_of_dif(double ax, double ay, double az, double bx, double by, double bz){
    return sqrt(pow_2(ax-bx) + pow_2(ay-by) + pow_2(az-bz));
}
double dot(double vx, double vy, double vz, double px, double py, double pz){
    //printf("%lf*%d + %lf*%d + %lf*%d\n", vx,px ,vy,py , vz,pz);
    return vx*px + vy*py + vz*pz;
}

double scale_I_KACh = 1.0, scale_I_pCa = 1.0, scale_I_bNa = 1.0, scale_I_bCa = 1.0, scale_I_NaCa = 1.0, scale_I_NaK = 1.0, scale_I_f = 1.0, scale_I_CaT = 1.0, scale_I_CaL = 1.0, scale_I_to = 1.0, scale_I_Kr = 1.0, scale_I_K1 = 1.0, scale_I_Kur = 1.0, scale_I_Ks = 1.0, scale_I_Na = 1.0;

// cell type defining variable
int cell_type = 1; // 1 = default
int SAN_cell = 0;
double IK1_v_shift = 0.0; // 0 = default

#pragma omp threadprivate(scale_I_KACh, scale_I_pCa, scale_I_bNa, scale_I_bCa, scale_I_NaCa, scale_I_NaK, scale_I_f, scale_I_CaT, scale_I_CaL, scale_I_to, scale_I_Kr, scale_I_K1, scale_I_Kur, scale_I_Ks, scale_I_Na, IK1_v_shift, cell_type, SAN_cell, C_m) 

// default (non AF) values for the AF dependent variables
int AF_model = 0;
double Ito_Act_Shift = 0.0;
double INa_Act_Shift = 0.0;
double fIRel = 1.0;
double GSR_leak = 1.0;
double fRyR = 0.0;
double RyRTauScale = 1.0;
double SERCA_Scale = 1.0;
void AF_model_selection(){// g_Cal, g_to,... same for all cells
    ////---------------------------------------------------------------------------------------------------------------////
    ////                                              Atrial Fibrillation                                              ////
    ////---------------------------------------------------------------------------------------------------------------////
    AF_model =  0; // 0 = no atrial fibrillation (default)
    // default (non AF) values for the AF dependent variables
    Ito_Act_Shift = 0.0;
    INa_Act_Shift = 0.0;
    fIRel = 1.0;
    GSR_leak = 1.0;
    fRyR = 0.0;
    RyRTauScale = 1.0;
    SERCA_Scale = 1.0;

    if (AF_model == 1) { // Bosch based model  // based on Henggui's paper. in cardiovascular research 2005.
		g_CaL *= 0.26;
		g_to *= 0.3;
		g_K1 *= 2.06;
		Ito_Act_Shift = 16;
		INa_Act_Shift = 1.6;
	}
	else if (AF_model == 2) { // Workman based model
		g_CaL *= 0.35;
		g_to *= 0.35;
		g_K1 *= 1.75;
	}
	else if (AF_model == 3) { // Wagoner based model
		g_CaL *= 0.37;
		g_to *= 0.34;
		g_Kur *= 0.51;
		g_K1 *= 2.06;
	}
	else if (AF_model == 4) { // Generic
		g_CaL *= 0.35;
		g_to *= 0.35;
		g_K1 *= 2.1;
		g_Kur *= 0.5;
        g_Ks *= 2;
		I_NaCamax *= 1.40;
        fIRel = 3.0;
		GSR_leak = 1.50; // Voigt et al. 2012
		fRyR = 2.5;
		RyRTauScale = 2.7;
		SERCA_Scale = 0.6;// voigt 2012 online figure III
		Ito_Act_Shift = 16;
		INa_Act_Shift = 1.6;
	}
    else if (AF_model == 5){ // Mechanisms of Human Atrial Fibrillation Initiation (2012), David E. Krumme,  Natalia A. Trayanova...
        g_CaL *= 0.3;
        g_to *= 0.5;
        g_Kur *= 0.5;
    }
}


void cell_type_selection(int sparse_index){// scale_I_Cal, scale_I_to,... changes depending on cell
    C_m = 0.0;
    scale_I_CaT = 1.0; // SAN only
    scale_I_CaL = 1.0;
    scale_I_to = 1.0;
    scale_I_Kr = 1.0;
    scale_I_K1 = 1.0; // atria only
    scale_I_f = 1.0; // SAN only
    scale_I_NaK = 1.0;
    scale_I_Kur = 1.0;
    scale_I_Ks = 1.0;
    scale_I_Na = 1.0;
    scale_I_NaCa = 1.0;
    scale_I_bCa = 1.0; // atria only
    scale_I_bNa = 1.0; // atria only
    scale_I_pCa = 1.0; // atria only
    scale_I_KACh = 1.0;
    IK1_v_shift = 0.0;
    SAN_cell = 0;
    D1_var = 0.0;
    D2_var = 0.0;

    D1_atria_var = 0.25;
    D2_atria_var = D1_atria_var*0.1;

    D1_SAN_var = 0.25;// temp 0.0012;
    D2_SAN_var = D1_SAN_var*0.1;

    if (cell_map_sparse[sparse_index]  == 1) {// normal atria, default Haibo 2017 model
        D1_var = D1_atria_var;
        D2_var = D2_atria_var;
        scale_I_Na *= 1.0*(1 - 0.2*0.0);
	}
    else if (cell_map_sparse[sparse_index]  == 2){ // generic fibrotic atrial tissue
        D1_var = D1_atria_var*0.125;
        D2_var = D2_atria_var*0.125;
    }
    else if (cell_map_sparse[sparse_index]  == 100){ // SAN wall
        D1_var = D2_atria_var*0.001;//D2_atria_var*0.05;//D2_atria_var*0.01; // D1 = D2 because wall should be isometric
        D2_var = D2_atria_var*0.001;//D2_atria_var*0.05;//D2_atria_var*0.01;
        scale_I_CaT *= 0.0; // SAN only
        scale_I_CaL *= 0.0;
        scale_I_to *= 0.0;
        scale_I_Kr *= 0.0;
        scale_I_K1 *= 0.0; // atria only
        scale_I_f *= 0.0; // SAN only
        scale_I_NaK *= 0.0;
        scale_I_Kur *= 0.0;
        scale_I_Ks *= 0.0;
        scale_I_Na *= 0.0;
        scale_I_NaCa *= 0.0;
        scale_I_bCa *= 0.0; // atria only
        scale_I_bNa *= 0.0; // atria only
        scale_I_pCa *= 0.0; // atria only
        scale_I_KACh *= 0.0;
    }
	
    else if (cell_map_sparse[sparse_index]  == 21){ // SAN head
        SAN_cell = 1;
        D1_var = D1_SAN_var*0.4*0.5;
        D2_var = D2_SAN_var*0.4*0.5;
        scale_I_K1 *= 0.0;
        scale_I_Na *= 2.0*(1 - 0.0);
        scale_I_f *= 0.5;
        scale_I_CaL *= 1.0;
        scale_I_KACh *= 1.0*0.1;
    }
    else if (cell_map_sparse[sparse_index]  == 22){ // SAN central
        SAN_cell = 1;
        D1_var = D1_SAN_var*0.4*0.5*0.7;
        D2_var = D2_SAN_var*0.4*0.5*0.7; // temp
        scale_I_K1 *= 0.0;

        scale_I_Na *= 1.0*(1 - 0.0);
        scale_I_f *= 1.0;
        scale_I_CaL *= 1.0;
        scale_I_KACh *= 1.0;
    }
    else if (cell_map_sparse[sparse_index]  == 23){ // SAN tail
        SAN_cell = 1;
        D1_var = D1_SAN_var*0.4*0.5;
        D2_var = D2_SAN_var*0.4*0.5;
        scale_I_K1 *= 0.0;
        scale_I_Na *= 2.0*(1 - 0.0);
        scale_I_f *= 0.5;
        scale_I_CaL *= 1.0;
        scale_I_KACh *= 1.0*0.1;
    }
    else if (cell_map_sparse[sparse_index]  == 24){ // SAN pathway
        SAN_cell = 1;
        scale_I_K1 *= 0.2;
        scale_I_Na *= 5.0*(1 - 0.0);
        scale_I_f *= 0.3;
        scale_I_CaL *= 1.0;
        scale_I_KACh *= 1.0;
        if (geo_y_index[sparse_index] <= 150){
            // left two pathways
            D1_var = D1_SAN_var*0.4; // temp
            D2_var = D2_SAN_var*0.4; // temp
        }
        else{ 
            // right three pathways
            D1_var = D1_SAN_var*0.4; // temp
            D2_var = D2_SAN_var*0.4; // temp
        } 
    }
    if(SAN_cell == 1){
        C_m = 57.0; // pF
    }
    else{
        C_m = 100.0;           // Membrane capacitance (pF)
    }
} 

#pragma omp threadprivate(D1_cube, D2_cube, fx_cube, fy_cube, fz_cube, V_cube)
#pragma omp threadprivate(I_Na, I_K1, I_to, I_Kur, I_Kr, I_Ks, I_CaL, I_CaT, I_f, I_NaK, I_NaCa, I_bCa, I_bNa, I_pCa, I_KACh, I_ion)
#pragma omp threadprivate(E_K, E_Na, E_Ca)
#pragma omp threadprivate(fx, fy, fz, dfxdx, dfxdy, dfxdz, dfydx, dfydy, dfydz, dfzdx, dfzdy, dfzdz)
#pragma omp threadprivate(D000, Dp00, Dm00, D0p0, D0m0, D00p, D00m, Dpp0, Dpm0, Dmp0, Dmm0, Dp0p, Dp0m, Dm0p, Dm0m, D0pp, D0pm, D0mp, D0mm)
#pragma omp threadprivate(alpha_m, beta_m, m_infinity, tau_m, alpha_h, beta_h, h_infinity, tau_h, alpha_j, beta_j, j_infinity, tau_j)
#pragma omp threadprivate(oa_infinity, tau_oa, oi_infinity, tau_oi)  
#pragma omp threadprivate(ua_infinity, tau_ua, ui_infinity, tau_ui) 
#pragma omp threadprivate(alpha_xr, beta_xr, tau_xr, xr_infinity, tau_paF, tau_paS, pa_infinity)
#pragma omp threadprivate(alpha_xs, beta_xs, tau_xs, xs_infinity)
#pragma omp threadprivate(tau_d, d_infinity, tau_f, f_infinity, tau_fCa, fCa_infinity, I_siCa, I_siNa, I_siK, alpha_d, beta_d, adVm, bdVm)
#pragma omp threadprivate(dT_infinity, fT_infinity, tau_dT, tau_fT)
#pragma omp threadprivate(f_NaK)
#pragma omp threadprivate(di, dodo, k12, k14, k21, k23, k32, k34, k41, k43, x1, x2, x3, x4)
#pragma omp threadprivate(tau_y, y_infinity, I_f_Na, I_f_K)
#pragma omp threadprivate(alpha_a, beta_a, tau_a, a_infinity)
#pragma omp threadprivate(BSS, B_i_NJ, B_sr_NJ, B_sr_SS)
#pragma omp threadprivate(I_serca_SR_NJ, I_serca_bulk_NJ, I_serca_SR_SS, I_serca_bulk_SS, I_SS_NJ)
#pragma omp threadprivate(RyR_serca_SS, RyR_a_SS_infinity, RyR_o_SS_infinity, RyR_c_SS_infinity, I_rel_SS)
#pragma omp threadprivate(RyR_serca_NJ, RyR_a_NJ_infinity, RyR_o_NJ_infinity, RyR_c_NJ_infinity, I_rel_NJ)
#pragma omp threadprivate(I_SRleak_NJ, I_SRleak_SS, I_Ca_NJ, I_Ca_SS, I_SRCa_NJ, I_SRCa_SS) 
#pragma omp threadprivate(dCa_SR_NJ, dCa_SR_SS)
#pragma omp threadprivate(j_SRCarel, diff, kCaSR, koSRCa, kiSRCa, delta_R, delta_O, delta_I, delta_RI, P_tot)
#pragma omp threadprivate(delta_fTC, delta_fTMC, delta_fTMM , delta_fCMi , delta_fCMs , delta_fCQ)
#pragma omp threadprivate(j_Ca_dif, ACh_block_P_up, j_up, delta_Ca_i, delta_Ca_sub, j_tr, delta_Ca_nsr, delta_Ca_jsr)

int main(int argc, char **argv)
{
    /*int num_threads = NUM_THREADS;
    if(NUM_CELLS < num_threads){
        num_threads = NUM_CELLS;
    }
    omp_set_num_threads(num_threads);*/
    time_t real_time, begin_timer = time(NULL), end_timer;
	int i, j, k, X_dim, Y_dim, Z_dim, adj_index;
    int sparse_index, number_of_cells;
	char c;
	FILE *in;
	FILE *out;
	char *str;
	double d2udx2, d2udy2, d2udz2, d2udxdy, d2udxdz, d2udydz, dudx, dudy, dudz, du, gz;
	
    double V_scale_high;
    double V_scale_low;
    double V_scale_gradient;
    double V_scale_constant;
    int V_print_val;
    int low_memory_output;
    int output_in_integers;
    int decimal_places_in_output_floats;

    // matlab auto saving variables
    char geometry_file_name[300];

    int save_to_mat_file = 1; // temp
    int delete_raw_files_after_saving_mat_file;
    int make_figure_images_from_mat_file;
    if (save_to_mat_file == 1){
        delete_raw_files_after_saving_mat_file = 1;
        make_figure_images_from_mat_file = 0;
    }

    low_memory_output = 1; // 0 = off, 1 = on. much lower memory (in sparse geometries), but not human readable
    
    if(low_memory_output == 0){

        output_in_integers = 0; // uses less memory than the default floats
        
        if(output_in_integers == 0){
            decimal_places_in_output_floats = 10;
        }
    }

    char *c_file_name;
    c_file_name = malloc(300*sizeof(char));
    sprintf(c_file_name,"%s", __FILE__);
    for(i = strlen(c_file_name)-1; i >= 0; i--){
        if(c_file_name[i] == '.') {
            c_file_name[i] = '\0';
            break;
        }
    }

    long int count = 0; // used to number the file names
	long int increment = 0; // increments the main time loop
    int increment_divisor = 400; //temp //120; // determines how often information gets printed to an output file
    double t = 0.0; // time (ms)
    int time_between_stimuli = (int)round(340); //260; // amount of time that should occur before next stimulus is applied (ms)
    int time_of_stimulation = 0; // the time at which the stimulation is applied (ms)
    int delt = 0.0;//50; // change in time between stimuli (ms)
    double stimulus_voltage = 20.0; // voltage of the stimulus being applied (mV)
    int stimulus_counter = 0; // keeps count of how many stimuli have occured
    int stimulus_to_continue;

    int pacing_location = (int)round(1.0);

    int initial_conditions_from_files = 0; // 0 = no (fresh run), 1 = yes (continued run)
    disable_burst_pacing = 1; // 0 = have a stimulaus, 1 = no stimulus
    
    // ############################## naming and creating main output directory ##############################
    int auto_naming = 1; // sets the name of the output directory to be the same as the c file that was executed
    
    char *main_output_directory_name;
	main_output_directory_name = malloc(300*sizeof(char));
    if(auto_naming){
	    sprintf(main_output_directory_name, "%s_output", c_file_name); 
    }
    else{
        char manually_set_name[300] = "example_manually_set_name";    // sets the name of the output directory, manually set, not recommended, as you may 
        sprintf(main_output_directory_name, "%s", manually_set_name); // accedentally overwrite the results of previous experiments, if your forget to change this
    }
    struct stat st = {0};
	if(stat(main_output_directory_name, &st) == -1) {mkdir(main_output_directory_name, 0777);}
    char saved_states_original_directory_name[300] = "saved_states_original"; // directory to save the saved state values in the original run
    char *saved_states_original_directory;
    saved_states_original_directory = malloc(300*sizeof(char));
    sprintf(saved_states_original_directory, "%s/%s", main_output_directory_name, saved_states_original_directory_name);
    char saved_states_rerun_directory_name[300] = "saved_states_rerun"; // directory to save the saved state values in subsequent runs
    char *saved_states_rerun_directory;
    saved_states_rerun_directory = malloc(300*sizeof(char));
    sprintf(saved_states_rerun_directory, "%s/%s", main_output_directory_name, saved_states_rerun_directory_name);
    char main_run_directory_name[300] = "results"; // directory to save the simulation results
    char *main_run_directory;
    main_run_directory = malloc(300*sizeof(char));
    sprintf(main_run_directory, "%s/%s", main_output_directory_name, main_run_directory_name);
    if(stat(saved_states_original_directory, &st) == -1) {mkdir(saved_states_original_directory, 0777);}
    if(stat(saved_states_rerun_directory, &st) == -1) {mkdir(saved_states_rerun_directory, 0777);}
    if(stat(main_run_directory, &st) == -1) {mkdir(main_run_directory, 0777);}

    char *stimulus_file;
    if(initial_conditions_from_files){

        stimulus_to_continue = 2;

        str = malloc(300*sizeof(char));
        sprintf(str, "cd %s; unzip stimulus_%03d.zip; cd -", saved_states_original_directory, stimulus_to_continue);
        int statuss;
        statuss = system(str);
        free(str);

        str = malloc(300*sizeof(char));
        sprintf(str, "%s/stimulus_%03d.dat", saved_states_original_directory, stimulus_to_continue);
        in = fopen(str, "r"); if(in == NULL){printf("\nfile_not_found: %s\n", str);}
        fscanf(in, "%d %lf %lf %d %d %d %ld %ld\n", &stimulus_counter, &t, &stimulus_voltage,  &time_of_stimulation, &time_between_stimuli, &delt, &increment, &count);
        stimulus_counter--;

        sprintf(str, "%s/L%d_continued_from_beat_%03d", main_run_directory, pacing_location, stimulus_to_continue);
        if (stat(str, &st) == -1) {
            mkdir(str, 0777);
        }
        free(str);
    }
    else{
        str = malloc(300*sizeof(char));
        sprintf(str, "%s/L%d", main_run_directory, pacing_location);
        if (stat(str, &st) == -1) {
            mkdir(str, 0777);
        }
        free(str);
    }

    /* Model Constant values (from Table 1 of the CRN 1998 paper) */
    
    g_Na = 7.8;            // Maximal I_Na conductance  (nS/pF)
    g_Na_SAN = (1000.0/57.0)*(0.0223);
    g_K1 = 0.09;           // Maximal I_K1 conductance  (nS/pF)
    g_K1_SAN = 0.09;
    g_to = 0.75471 * 0.1962;         // Maximal I_to conductance  (nS/pF)
    g_to_SAN = (1000.0/57.0)*(0.0035);
    g_Kur = 0.9 * 0.006398;
    g_Kur_SAN = (1000.0/57.0)*(0.1539e-3);
    g_Kr = 0.8 * 0.029411765;         // Maximal I_Kr conductance  (nS/pF)
    g_Kr_SAN = (1000.0/57.0)*(0.00424);
    g_Ks = 0.8 * 0.12941176;          // Maximal I_Ks conductance  (nS/pF)
    g_Ks_SAN = (1000.0/57.0)*(0.00065);
    g_CaL = 0.75 * 0.1294;        // Maximal I_CaL conductance (nS/pF)
    P_CaL = (1000.0/57.0)*(0.4578);
    g_bCa = 1.4 * 0.001131;       // Maximal I_bCa conductance (nS/pF)
    g_bNa = 0.8 * 0.0006744375;      // Maximal I_bNa conductance (nS/pF)
    I_NaKmax = 1.4 * 0.59933874;       // Maximal I_NaK  (pA/pF)
    I_NaKmax_SAN = (1000.0/57.0)*(0.08105);
    sigma = (exp(Na_o/67.3)-1)/7.0;    // Na_o dependance parameter for f_NaK (dimensionless)
    I_NaCamax = 1600.0;    // Maximal I_NaCa (pA/pF)
    I_pCamax = 0.275;      // Maximal I_pCa  (pA/pF)
    ACh = 0.0;//1e-6;  //0.0; // temp        // Acetylcholine concentration (mM)
    E_K_SAN = (R*T/F)*log(K_o/K_i_SAN);
    E_Na_SAN = (R*T/F)*log(Na_o/Na_i_SAN);;
    E_mh_SAN = (R*T/F)*log((Na_o + 0.12*K_o)/(Na_i_SAN + 0.12*K_i_SAN));
    E_Ks_SAN = (R*T/F)*log((K_o + 0.12*Na_o)/(K_i_SAN + 0.12*Na_i_SAN));
    
    G_f = g_f/(K_o/(K_o+Km_f));
    G_f_K = G_f/(alpha+1.0);
    G_f_Na = alpha*G_f_K;
    g_f_Na = (1000.0/57.0)*G_f_Na*K_o/(K_o+Km_f); // (nS/pF)
    g_f_K = (1000.0/57.0)*G_f_K*K_o/(K_o+Km_f); // (nS/pF)
    g_KACh_SAN = (1000.0/57.0)*(0.00345);
    if(fabs(ACh) > 1e-9){
        ACh_shift_I_f = -1.0 - 9.898*pow(1.0*ACh, 0.618)/(pow(1.0*ACh, 0.618)+0.00122423);
        ACh_block_I_CaL = 0.31*ACh/(ACh+0.00009);
    }
    else{
        ACh_shift_I_f = 0.0;
        ACh_block_I_CaL = 0.0;
    }
    P_CaT = (1000.0/57.0)*0.04132;
    V_sub = 0.000000001*2.0*3.14159265358979*L_sub*(R_cell-L_sub/2.0)*L_cell;
    V_cell = 0.000000001*3.14159265358979*pow_2(R_cell)*L_cell;
    V_nsr = V_nsr_part*V_cell;
    V_i = V_i_part*V_cell-V_sub;
    V_jsr = V_jsr_part*V_cell;
    I_NaCamax_SAN = (1000.0/57.0)*3.343;
////---------------------------------------------------------------------------------------------------------------////
////                                              Atrial Fibrillation                                              ////
////---------------------------------------------------------------------------------------------------------------////
               AF_model_selection();// needs to be placed after cell paraneters are given their values                                
////---------------------------------------------------------------------------------------------------------------////

    int tm1, tm2, tm3, tm4, tm5, tm6, tm7, tm8, tm9, tm10, tm11, tm12, tm13, tm14, tm15, tm16, tm17, tm18; // unused variables
    int cell_map_num;
    sprintf(geometry_file_name, "SAN_1_raised_cen_hist_fibres.dat");
    in = fopen(geometry_file_name, "r"); if(in == NULL){printf("\ngeometry_file_not_found\n");}
    fscanf(in, "%d %d %d %d %d %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", &number_of_cells, &X_dim, &Y_dim, &Z_dim, &cell_map_num, &V_scale_high, &V_scale_low, &gz, &tm1, &tm2, &tm3, &tm4, &tm5, &tm6, &tm7, &tm8, &tm9, &tm10, &tm11, &tm12, &tm13, &tm14, &tm15, &tm16, &tm17, &tm18);
    for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
        fscanf(in, "%d %hd %hd %hd %hd %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", 
            &tm1, &geo_x_index[sparse_index], &geo_y_index[sparse_index], &geo_z_index[sparse_index], &cell_map_sparse[sparse_index],
            &fx_sparse[sparse_index], &fy_sparse[sparse_index], &fz_sparse[sparse_index],
                                        //  x,  y,  z 
            &geo_adj[sparse_index][0],  // -1, -1,  0
            &geo_adj[sparse_index][1],  // -1,  0, -1
            &geo_adj[sparse_index][2],  // -1,  0,  0
            &geo_adj[sparse_index][3],  // -1,  0,  1
            &geo_adj[sparse_index][4],  // -1,  1,  0
            &geo_adj[sparse_index][5],  //  0, -1, -1
            &geo_adj[sparse_index][6],  //  0, -1,  0
            &geo_adj[sparse_index][7],  //  0, -1,  1
            &geo_adj[sparse_index][8],  //  0,  0, -1
            &geo_adj[sparse_index][9],  //  0,  0,  1
            &geo_adj[sparse_index][10], //  0,  1, -1
            &geo_adj[sparse_index][11], //  0,  1,  0
            &geo_adj[sparse_index][12], //  0,  1,  1
            &geo_adj[sparse_index][13], //  1, -1,  0
            &geo_adj[sparse_index][14], //  1,  0, -1
            &geo_adj[sparse_index][15], //  1,  0,  0
            &geo_adj[sparse_index][16], //  1,  0,  1
            &geo_adj[sparse_index][17]);//  1,  1,  0         
    }
    fclose(in);
    
    V_scale_gradient = 127.0/(V_scale_high - V_scale_low);
    V_scale_constant = -V_scale_low*V_scale_gradient;

    //*
    for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
        cell_type_selection(sparse_index);
        D1m_sparse[sparse_index] = D1_var;
        D2m_sparse[sparse_index] = D2_var;
    }
    //*/

    //for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
    //    printf("cell_num = %d: D1=%lf, D2=%lf\n", sparse_index, D1m_sparse[sparse_index], D2m_sparse[sparse_index]);
    //}

    //setting initial conditions for all points
    if(initial_conditions_from_files){
        double temp_1, temp_4;
        long int temp_2, temp_3;
        int nothing_1, nothing_2, nothing_3, nothing_4;
        str = malloc(300*sizeof(char));
        sprintf(str, "%s/stimulus_%03d.dat", saved_states_original_directory, stimulus_to_continue);
        in = fopen(str, "r"); if(in == NULL){printf("\nfile_not_found: %s\n", str);}
        fscanf(in, "%d %lf %lf %d %d %d %ld %ld\n", &nothing_1, &temp_1, &temp_4, &nothing_2, &nothing_3, &nothing_4, &temp_2, &temp_3); // properly read in at line ~626
        for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
            fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                &V_sparse[sparse_index], &Na_i_sparse[sparse_index], &K_i_sparse[sparse_index], &Ca_i_SS_sparse[sparse_index], &Ca_sr_SS_sparse[sparse_index], &Ca_serca_SS_sparse[sparse_index], &m_sparse[sparse_index], 
                &h_sparse[sparse_index], &j_sparse[sparse_index], &oa_sparse[sparse_index], &oi_sparse[sparse_index], &ua_sparse[sparse_index], &ui_sparse[sparse_index], &xr_sparse[sparse_index], &xs_sparse[sparse_index],
                &d_sparse[sparse_index], &f_sparse[sparse_index], &fCa_sparse[sparse_index], &Ca_i_NJ_sparse[sparse_index], &Ca_sr_NJ_sparse[sparse_index], &Ca_serca_NJ_sparse[sparse_index],
                &RyR_o_SS_sparse[sparse_index], &RyR_c_SS_sparse[sparse_index], &RyR_a_SS_sparse[sparse_index], &RyR_o_NJ_sparse[sparse_index], &RyR_c_NJ_sparse[sparse_index], &RyR_a_NJ_sparse[sparse_index]);
        }
        fclose(in);
        sprintf(str, "cd %s; rm stimulus_%03d.dat; cd -", saved_states_original_directory, stimulus_to_continue);
        int statuss;
        statuss = system(str);
        free(str);
    }
    else{
        

        int forced_reentrant = 0;
        if(forced_reentrant){
            disable_burst_pacing = 1;
            double rmp_plane[5], stim_plane[5], refr_plane[5]; // ax + by + cz + d = 0 (5th index is for type of equality to 0)
            double p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, pa, pb, pc, pd, pe, dif_norm;
            double rest_p_v[3], stim_p_v[3], refr_p_v[3]; // normal vectors for planes
            int s1_s2_points[20][3][3]; // 20 locations, and 3 points with 3 coordinates each
            char IC_string_rmp[500], IC_string_stim[500], IC_string_refr[500], IC_string_to_apply[500];
            char *line_buf = NULL;
            size_t line_buf_size = 0;
            int line_count = 0;
            ssize_t line_size;
            in = fopen("s1s2.dat", "r"); if(in == NULL){printf("\ns1 s2 set up file_not_found\n");}
            line_size = getline(&line_buf, &line_buf_size, in);
            sprintf(IC_string_rmp, "%s", line_buf);
            line_size = getline(&line_buf, &line_buf_size, in);
            sprintf(IC_string_stim, "%s", line_buf);
            line_size = getline(&line_buf, &line_buf_size, in);
            sprintf(IC_string_refr, "%s", line_buf);
            while (line_size >= 0){
                line_size = getline(&line_buf, &line_buf_size, in);
                if(line_size > 0){
                    sscanf(line_buf, "%d %d %d %d %d %d %d %d %d\n", 
                        &s1_s2_points[line_count][1][1], &s1_s2_points[line_count][1][2], &s1_s2_points[line_count][1][3],
                        &s1_s2_points[line_count][2][1], &s1_s2_points[line_count][2][2], &s1_s2_points[line_count][2][3],
                        &s1_s2_points[line_count][3][1], &s1_s2_points[line_count][3][2], &s1_s2_points[line_count][3][3]
                    );
                    printf("%s",line_buf);
                }
                line_count++;
            }
            fclose(in);

            int s1_s2_location = 1;
            p1x=1.0*s1_s2_points[s1_s2_location-1][1][1]; p1y=1.0*s1_s2_points[s1_s2_location-1][1][2]; p1z=1.0*s1_s2_points[s1_s2_location-1][1][3];
            p2x=1.0*s1_s2_points[s1_s2_location-1][2][1]; p2y=1.0*s1_s2_points[s1_s2_location-1][2][2]; p2z=1.0*s1_s2_points[s1_s2_location-1][2][3];
            p3x=1.0*s1_s2_points[s1_s2_location-1][3][1]; p3y=1.0*s1_s2_points[s1_s2_location-1][3][2]; p3z=1.0*s1_s2_points[s1_s2_location-1][3][3];
            
            // define IC planes
            // refractory plane
            //printf("%lf %lf %lf %lf %lf %lf\n", p3x, p3y, p3z, p2x, p2y, p2z);
            dif_norm = norm_of_dif(p3x, p3y, p3z, p2x, p2y, p2z);
            pa = (p3x-p2x)/dif_norm; pb = (p3y-p2y)/dif_norm; pc = (p3z-p2z)/dif_norm;
            pd = -1.0*dot(pa, pb, pc, p3x, p3y, p3z);
            pe = 1; if((pa*p2x + pb*p2y + pc*p2z + pd)*pe >= 0.0){pe = -1;}
            refr_plane[0] = pa; refr_plane[1] = pb; refr_plane[2] = pc; refr_plane[3] = pd; refr_plane[4] = pe;

            // stim plane
            dif_norm = norm_of_dif(p2x, p2y, p2z, p1x, p1y, p1z);
            pa = (p2x-p1x)/dif_norm; pb = (p2y-p1y)/dif_norm; pc = (p2z-p1z)/dif_norm;
            pd = -1.0*dot(pa, pb, pc, p2x, p2y, p2z);
            pe = 1; if((pa*p1x + pb*p1y + pc*p1z + pd)*pe >= 0.0){pe = -1;}
            stim_plane[0] = pa; stim_plane[1] = pb; stim_plane[2] = pc; stim_plane[3] = pd; stim_plane[4] = pe;

            // resting plane
            dif_norm = norm_of_dif(p3x, p3y, p3z, p2x, p2y, p2z);
            pa = (p3x-p2x)/dif_norm; pb = (p3y-p2y)/dif_norm; pc = (p3z-p2z)/dif_norm;
            pd = -1.0*dot(pa, pb, pc, p2x, p2y, p2z);
            pe = 1; if((pa*p3x + pb*p3y + pc*p3z + pd)*pe >= 0.0){pe = -1;}
            rmp_plane[0] = pa; rmp_plane[1] = pb; rmp_plane[2] = pc; rmp_plane[3] = pd; rmp_plane[4] = pe;

            // make IC mask
            for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                if(!((refr_plane[0]*geo_x_index[sparse_index] + refr_plane[1]*geo_y_index[sparse_index] + refr_plane[2]*geo_z_index[sparse_index] + refr_plane[3])*refr_plane[4] >= 0.0)){
                    if(!((rmp_plane[0]*geo_x_index[sparse_index] + rmp_plane[1]*geo_y_index[sparse_index] + rmp_plane[2]*geo_z_index[sparse_index] + rmp_plane[3])*rmp_plane[4] >= 0.0)
                    && !((stim_plane[0]*geo_x_index[sparse_index] + stim_plane[1]*geo_y_index[sparse_index] + stim_plane[2]*geo_z_index[sparse_index] + stim_plane[3])*stim_plane[4]  >= 0.0)   ){
                        IC_map[sparse_index] = 2; // stimulus
                    }
                    else{
                        IC_map[sparse_index] = 1; // resting membrane potential
                    }
                }
                else{
                    IC_map[sparse_index] = 3; // refractory
                }
            }

            // exclude regions unconnected to point_3 in an s1s2 trio

            // assign ICs according to IC mask labels
            for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                if     (IC_map[sparse_index] == 1){sprintf(IC_string_to_apply, "%s", IC_string_rmp);}
                else if(IC_map[sparse_index] == 2){sprintf(IC_string_to_apply, "%s", IC_string_stim);}
                else if(IC_map[sparse_index] == 3){sprintf(IC_string_to_apply, "%s", IC_string_refr);}
                sscanf(IC_string_to_apply, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                    &V_sparse[sparse_index], &Na_i_sparse[sparse_index], &K_i_sparse[sparse_index], &Ca_i_SS_sparse[sparse_index], &Ca_sr_SS_sparse[sparse_index], &Ca_serca_SS_sparse[sparse_index], &m_sparse[sparse_index], 
                    &h_sparse[sparse_index], &j_sparse[sparse_index], &oa_sparse[sparse_index], &oi_sparse[sparse_index], &ua_sparse[sparse_index], &ui_sparse[sparse_index], &xr_sparse[sparse_index], &xs_sparse[sparse_index],
                    &d_sparse[sparse_index], &f_sparse[sparse_index], &fCa_sparse[sparse_index], &Ca_i_NJ_sparse[sparse_index], &Ca_sr_NJ_sparse[sparse_index], &Ca_serca_NJ_sparse[sparse_index],
                    &RyR_o_SS_sparse[sparse_index], &RyR_c_SS_sparse[sparse_index], &RyR_a_SS_sparse[sparse_index], &RyR_o_NJ_sparse[sparse_index], &RyR_c_NJ_sparse[sparse_index], &RyR_a_NJ_sparse[sparse_index]);
            }
        }
        else{ // default ICs
            for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                if((cell_map_sparse[sparse_index] >= 21) && (cell_map_sparse[sparse_index] <= 29)){ // SAN cells
                    V_sparse[sparse_index] = -58.7607;
                    m_sparse[sparse_index] = 0.1181030;
                    h_sparse[sparse_index] = 0.0425886;
                    oa_sparse[sparse_index] = 0.4592; // q from Fabbri
                    oi_sparse[sparse_index] = 0.00550695; // r from Fabbri

                    ua_sparse[sparse_index] = 0.00401424; // r_Kur from Fabbri
                    ui_sparse[sparse_index] = 0.835394; // s_Kur from Fabbri
                    paF_sparse[sparse_index] = 0.0150839;
                    paS_sparse[sparse_index] = 0.66063;
                    xr_sparse[sparse_index] = 0.853491; // piy from Fabbri

                    xs_sparse[sparse_index] = 0.167373; // n from Fabbri

                    d_sparse[sparse_index] = 5.83515e-05;
                    f_sparse[sparse_index] = 0.513361;
                    fCa_sparse[sparse_index] = 0.636286;

                    dT_sparse[sparse_index] = 0.0236757;
                    fT_sparse[sparse_index] = 0.383249;

                    y_sparse[sparse_index] = 0.00295033;
                    a_sparse[sparse_index] = 0.00204474; 
  
                    RyR_o_SS_sparse[sparse_index]     = 1.05007e-09;      // I (dimensionless) (in Ca_SR_release)
                    RyR_o_NJ_sparse[sparse_index]     = 3.94166e-09;       // O (dimensionless) (in Ca_SR_release)
                    RyR_a_SS_sparse[sparse_index]     = 0.7896370;         // R1 (dimensionless) (R in Ca_SR_release)
                    RyR_a_NJ_sparse[sparse_index]     = 0.210362;         // RI (dimensionless) (in Ca_SR_release)
                    RyR_c_SS_sparse[sparse_index]     = 0.297126;        // fCMi (dimensionless) (in Ca_buffering)
                    RyR_c_NJ_sparse[sparse_index]     = 0.158989;        // fCMs (dimensionless) (in Ca_buffering)
                    fCQ_sparse[sparse_index]          = 0.0978994;         // fCQ (dimensionless) (in Ca_buffering)
                    fTC_sparse[sparse_index]          = 0.0270936;         // fTC (dimensionless) (in Ca_buffering)
                    Ca_serca_SS_sparse[sparse_index]  = 0.346474;        // fTMC (dimensionless) (in Ca_buffering)
                    Ca_serca_NJ_sparse[sparse_index]  = 0.577331;        // fTMM (dimensionless) (in Ca_buffering)
                    Ca_sr_SS_sparse[sparse_index]     = 0.27698;          // Ca_jsr (millimolar) (in Ca_dynamics)
                    Ca_sr_NJ_sparse[sparse_index]     = 0.419885;          // Ca_nsr (millimolar) (in Ca_dynamics)
                    Ca_i_SS_sparse[sparse_index]      = 6.2441e-05;       // Ca_sub (millimolar) (in Ca_dynamics)
                    Ca_i_NJ_sparse[sparse_index]      = 1.39073e-04;           // Cai (millimolar) (in Ca_dynamics)

                }
                else{
                    V_sparse[sparse_index] = -77.13255836;              // membrane voltage
                    m_sparse[sparse_index] = 0.005631819916;
                    h_sparse[sparse_index] = 0.9168420281;
                    j_sparse[sparse_index] = 0.9380183564;
                    d_sparse[sparse_index] = 0.0002267161277;           // I_CaL gates
                    f_sparse[sparse_index] = 0.9354212881;           
                    fCa_sparse[sparse_index] = 0.7270823666;         
                    xr_sparse[sparse_index] = 0.00157418133;            // I_Kr gate
                    xs_sparse[sparse_index] = 0.02225979641;            // I_Ks gates
                    oa_sparse[sparse_index] = 0.01223706011;            // I_to gates
                    oi_sparse[sparse_index] = 0.8849139842;
                    ua_sparse[sparse_index] = 0.0002417881801;          // I_Kur gates
                    ui_sparse[sparse_index] = 0.9517278864;           
                    Na_i_sparse[sparse_index] = 10.30397012;            // Na_i_sparse[sparse_index]
                    K_i_sparse[sparse_index] = 131.867138;              // ki
                    Ca_i_NJ_sparse[sparse_index] = 0.0001403133065;     // Calcium handling
                    Ca_i_SS_sparse[sparse_index] = 0.0001313595105;
                    Ca_sr_NJ_sparse[sparse_index] = 0.9892411621;
                    Ca_sr_SS_sparse[sparse_index] = 0.9779168037;
                    Ca_serca_NJ_sparse[sparse_index] = 0.00958584702;
                    Ca_serca_SS_sparse[sparse_index] = 0.009386941118;
                    RyR_o_SS_sparse[sparse_index] = 0.000456694441;
                    RyR_c_SS_sparse[sparse_index] = 0.9571976502;
                    RyR_a_SS_sparse[sparse_index] = 0.1419271573;
                    RyR_o_NJ_sparse[sparse_index] = 0.0003667999273;
                    RyR_c_NJ_sparse[sparse_index] = 0.9790238698;
                    RyR_a_NJ_sparse[sparse_index] = 0.2977219443;
                }
            }
        }

    }
    for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
        V_new_sparse[sparse_index] = V_sparse[sparse_index];
    }

    //V_sparse[62] = 20.0; // temp

    //while (time_between_stimuli >= 210) {
    while (t < 3000.0) { //temp
    //while (t < 1.1) { //temp
        if (fabs(time_of_stimulation - t) <= 1.1*dt){
            /*
            stimulus_counter++;
            stimulus_file = malloc(300*sizeof(char));
            printf("stimulus_%03d\n", stimulus_counter);
            if(initial_conditions_from_files == 1){
                sprintf(stimulus_file, "%s/stimulus_%03d.dat", saved_states_rerun_directory, stimulus_counter);
            }
            else{
                sprintf(stimulus_file, "%s/stimulus_%03d.dat", saved_states_original_directory, stimulus_counter);
            }
            out = fopen(stimulus_file, "wt");
            fprintf(out, "%d %lf %lf %d %d %d %ld %ld\n", stimulus_counter, t, stimulus_voltage,  time_of_stimulation, time_between_stimuli, delt, increment, count);
            for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                fprintf(out,"%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
                    V_sparse[sparse_index], Na_i_sparse[sparse_index], K_i_sparse[sparse_index], Ca_i_SS_sparse[sparse_index], Ca_sr_SS_sparse[sparse_index], Ca_serca_SS_sparse[sparse_index], m_sparse[sparse_index], h_sparse[sparse_index],
                    j_sparse[sparse_index], oa_sparse[sparse_index], oi_sparse[sparse_index], ua_sparse[sparse_index], ui_sparse[sparse_index], xr_sparse[sparse_index], xs_sparse[sparse_index], d_sparse[sparse_index], f_sparse[sparse_index], fCa_sparse[sparse_index],
                    Ca_i_NJ_sparse[sparse_index], Ca_sr_NJ_sparse[sparse_index], Ca_serca_NJ_sparse[sparse_index], RyR_o_SS_sparse[sparse_index], RyR_c_SS_sparse[sparse_index],
                    RyR_a_SS_sparse[sparse_index], RyR_o_NJ_sparse[sparse_index], RyR_c_NJ_sparse[sparse_index], RyR_a_NJ_sparse[sparse_index]);
            }      
            fclose(out);
            str = malloc(300*sizeof(char));
            if(initial_conditions_from_files){
                sprintf(str, "cd %s; zip stimulus_%03d.zip stimulus_%03d.dat; rm stimulus_%03d.dat; cd -", saved_states_rerun_directory, stimulus_counter, stimulus_counter, stimulus_counter);
            }
            else{
                sprintf(str, "cd %s; zip stimulus_%03d.zip stimulus_%03d.dat; rm stimulus_%03d.dat; cd -", saved_states_original_directory, stimulus_counter, stimulus_counter, stimulus_counter);
            }
            int statuss;
            statuss = system(str);
            free (str);*/

            
            // pacing_location_x is the Y coordinate from matlab graphs, and pacing_location_y is the X coordinate
            int pacing_location_x, pacing_location_y, pacing_location_z;
            if     (pacing_location ==  1) {pacing_location_x =  91;    pacing_location_y =  212;    pacing_location_z =  113;}
            else if(pacing_location ==  2) {pacing_location_x =  73;    pacing_location_y =  212;    pacing_location_z =   80;}
            else if(pacing_location ==  3) {pacing_location_x =  45;    pacing_location_y =  212;    pacing_location_z =   55;}
            else if(pacing_location ==  4) {pacing_location_x =  80;    pacing_location_y =  156;    pacing_location_z =  121;}
            else if(pacing_location ==  5) {pacing_location_x =  54;    pacing_location_y =  156;    pacing_location_z =   93;}
            else if(pacing_location ==  6) {pacing_location_x =  31;    pacing_location_y =  156;    pacing_location_z =   63;}
            else if(pacing_location ==  7) {pacing_location_x =  93;    pacing_location_y =   93;    pacing_location_z =  124;}
            else if(pacing_location ==  8) {pacing_location_x =  62;    pacing_location_y =   93;    pacing_location_z =   99;}
            else if(pacing_location ==  9) {pacing_location_x =  27;    pacing_location_y =   93;    pacing_location_z =   73;}
            else if(pacing_location == 10) {pacing_location_x =  60;    pacing_location_y =  198;    pacing_location_z =  71;}
            else if(pacing_location == 11) {pacing_location_x =  53;    pacing_location_y =  191;    pacing_location_z =  70;}
            else if(pacing_location == 12) {pacing_location_x =  67;    pacing_location_y =  123;    pacing_location_z =  115;}
            else if(pacing_location == 13) {pacing_location_x =  45;    pacing_location_y =  123;    pacing_location_z =  79;}
            else{printf("\ninvalid pacing location.\n");}

            if(disable_burst_pacing == 0){
                /*
                for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                    if(abs(geo_x_index[sparse_index] - pacing_location_x)  <= 3){
                        if(abs(geo_y_index[sparse_index] - pacing_location_y)  <= 3){
                            if(abs(geo_z_index[sparse_index] - pacing_location_z)  <= 3){
                                V_sparse[sparse_index] = stimulus_voltage;
                            }
                        }
                    }
                }
                //*/
                //*
                int xi = 264-25, xf = 264+25;
                int yi = 136-25, yf = 136+25; 
                int zi = 2, zf = 2;
                for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                    if((geo_x_index[sparse_index] >= xi) && (geo_x_index[sparse_index] <= xf)){
                        if((geo_y_index[sparse_index] >= yi) && (geo_y_index[sparse_index] <= yf)){
                            if((geo_z_index[sparse_index] >= zi) && (geo_z_index[sparse_index] <= zf)){
                                V_sparse[sparse_index] = stimulus_voltage;
                            }
                        }
                    }
                }
                //*/
            }

            //* set the time intervals to apply the stimulus after a few timesteps after the first stimulus is applied
            if (time_between_stimuli <= 500){
                delt = 10;
            }
            if (time_between_stimuli <= 300){
                delt = 10;
            }
            if (time_between_stimuli <= 200){
                delt = 10;
            }
            if (time_between_stimuli <= 100){
                delt = 5;
            }
            time_between_stimuli -=  delt;
            //*/
            time_of_stimulation += time_between_stimuli; 
	    }

        if ((increment % increment_divisor == 0)) {
            str = malloc (300*sizeof(char));
            if(initial_conditions_from_files){
                sprintf(str, "%s/L%d_continued_from_beat_%03d/rerun_%05ld.vtk", main_run_directory, pacing_location, stimulus_to_continue, count);
            }
            else{
                sprintf(str, "%s/L%d/a_%05ld.vtk",main_run_directory, pacing_location, count); //output
            }
            time(&real_time);
            printf("%s, File %ld, t=%.3lf ms, PCL:%d, Date & Time: %s", c_file_name, count, t, time_between_stimuli, ctime(&real_time));
            count++;
            out = fopen (str, "wt");
            if(low_memory_output == 1){
                for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
                    V_print_val = (int)round(V_sparse[sparse_index]*V_scale_gradient + V_scale_constant);
                    if(V_print_val > 127){
                        V_print_val = 127;
                    }
                    else if(V_print_val < 0){
                        V_print_val = 0;
                    }
                    fprintf(out, "%c", V_print_val);
                }
            }
            else{
                for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){ // need to modify matlab reader file
                    if(output_in_integers){
                        fprintf (out, "%d\n", (int)round(V_sparse[sparse_index]));
                    }
                    else{
                        fprintf (out, "%.*lf\n", decimal_places_in_output_floats, V_sparse[sparse_index]);
                    }
                }
            }
            fclose (out);
            free (str);
        }


        //V_sparse[62] = 10.0; // temp

        #pragma omp parallel for private(i, j, k, adj_index, sparse_index, d2udx2, d2udy2, d2udz2, d2udxdy, d2udxdz, d2udydz, dudx, dudy, dudz, du) 
        for(sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
            if(cell_map_sparse[sparse_index]  == 100){
                V_new_sparse[sparse_index] = -62.5; // forced voltage at walls
            }
            else{
            cell_type_selection(sparse_index);
            for(i = 0; i <= 2; i++){
                for(j = 0; j <= 2; j++){
                    for(k = 0; k <= 2; k++){
                        V_cube[i][j][k]=0.0; D1_cube[i][j][k]=0.0; D2_cube[i][j][k]=0.0; fx_cube[i][j][k]=0.0; fy_cube[i][j][k]=0.0; fz_cube[i][j][k]=0.0;
                    }
                }
            }
            adj_index = 0;
            for(i = 0; i <= 2; i++){
                for(j = 0; j <= 2; j++){
                    for(k = 0; k <= 2; k++){
                        if(((i-1)*(j-1)*(k-1) == 0) && !((i==1)&&(j==1)&&(k==1))){ // 18 closest neighbours
                            if(geo_adj[sparse_index][adj_index] != -1){ // neighbour cell
                                D1_cube[i][j][k] = D1m_sparse[ geo_adj[sparse_index][adj_index] ]; D2_cube[i][j][k] = D2m_sparse[ geo_adj[sparse_index][adj_index] ];
                                fx_cube[i][j][k] =  fx_sparse[ geo_adj[sparse_index][adj_index] ]; fy_cube[i][j][k] =  fy_sparse[ geo_adj[sparse_index][adj_index] ];
                                V_cube[i][j][k]  =   V_sparse[ geo_adj[sparse_index][adj_index] ]; fz_cube[i][j][k] =  fz_sparse[ geo_adj[sparse_index][adj_index] ];
                            }
                            adj_index++;
                        }
                        else{ 
                            if((i==1)&&(j==1)&&(k==1)){  // current cell
                                fx_cube[i][j][k] =  fx_sparse[sparse_index]; fy_cube[i][j][k] =  fy_sparse[sparse_index]; fz_cube[i][j][k] =  fz_sparse[sparse_index];
                                V_cube[i][j][k]  =   V_sparse[sparse_index]; D1_cube[i][j][k] = D1m_sparse[sparse_index]; D2_cube[i][j][k] = D2m_sparse[sparse_index];
                            }
                        }
                    }
                }
            }

            dfxdx = fx_cube[2][1][1] - fx_cube[0][1][1]; // dfx/dx,( 1/(2*dx) is accounted for in final calculation)
            dfxdy = fx_cube[1][2][1] - fx_cube[1][0][1]; // dfx/dy
            dfxdz = fx_cube[1][1][2] - fx_cube[1][1][0]; // dfx/dz
            dfydx = fy_cube[2][1][1] - fy_cube[0][1][1]; // dfy/dx
            dfydy = fy_cube[1][2][1] - fy_cube[1][0][1]; // dfy/dy
            dfydz = fy_cube[1][1][2] - fy_cube[1][1][0]; // dfy/dz
            dfzdx = fz_cube[2][1][1] - fz_cube[0][1][1]; // dfz/dx
            dfzdy = fz_cube[1][2][1] - fz_cube[1][0][1]; // dfz/dy
            dfzdz = fz_cube[1][1][2] - fz_cube[1][1][0]; // dfz/dz
            
            D000 = D1_cube[1][1][1];
            Dp00 = D1_cube[2][1][1]; D0p0 = D1_cube[1][2][1]; D00p = D1_cube[1][1][2];
            Dm00 = D1_cube[0][1][1]; D0m0 = D1_cube[1][0][1]; D00m = D1_cube[1][1][0];
            Dpp0 = D1_cube[2][2][1]; Dp0p = D1_cube[2][1][2]; D0pp = D1_cube[1][2][2];
            Dmm0 = D1_cube[0][0][1]; Dm0m = D1_cube[0][1][0]; D0mm = D1_cube[1][0][0];
            Dpm0 = D1_cube[2][0][1]; Dp0m = D1_cube[2][1][0]; D0pm = D1_cube[1][2][0];
            Dmp0 = D1_cube[0][2][1]; Dm0p = D1_cube[0][1][2]; D0mp = D1_cube[1][0][2];
            fx = fx_cube[1][1][1]; fy = fy_cube[1][1][1]; fz = fz_cube[1][1][1];
            //fx = 0.0; fy = 0.0; fz = 0.0; // temp

            du=0.0; // temp
            //*
            min_of_2(D000,D2_cube[1][1][1],  Dp00,D2_cube[2][1][1]);
            du += (V_cube[2][1][1] - V_cube[1][1][1])*((D1-D2)*fx*fx + D2)*4;
            //if((geo_x_index[sparse_index]==319)&&(geo_y_index[sparse_index]==149)&&(geo_z_index[sparse_index]==2)){
            //    printf("du=%lf %lf %lf %lf %lf %lf\n", (V_cube[2][1][1] - V_cube[1][1][1])*((D1-D2)*fx*fx + D2)*4, V_cube[2][1][1], V_cube[1][1][1], D1, D2, fx);
                //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", V_cube[0][0][1], V_cube[0][1][1], V_cube[0][2][1], V_cube[1][0][1], V_cube[1][1][1], V_cube[1][2][1], V_cube[2][0][1], V_cube[2][1][1], V_cube[2][2][1]);
            //}
            min_of_2(D000,D2_cube[1][1][1],  Dm00,D2_cube[0][1][1]);
            du += (V_cube[0][1][1] - V_cube[1][1][1])*((D1-D2)*fx*fx + D2)*4;
            min_of_2(D000,D2_cube[1][1][1],  D0p0,D2_cube[1][2][1]);
            du += (V_cube[1][2][1] - V_cube[1][1][1])*((D1-D2)*fy*fy + D2)*4;
            min_of_2(D000,D2_cube[1][1][1],  D0m0,D2_cube[1][0][1]);
            du += (V_cube[1][0][1] - V_cube[1][1][1])*((D1-D2)*fy*fy + D2)*4;
            //min_of_2(D000,D2_cube[1][1][1],  D00p,D2_cube[1][1][2]);
            //du += (V_cube[1][1][2] - V_cube[1][1][1])*((D1-D2)*fz*fz + D2)*4;
            //min_of_2(D000,D2_cube[1][1][1],  D00m,D2_cube[1][1][0]);
            //du += (V_cube[1][1][0] - V_cube[1][1][1])*((D1-D2)*fz*fz + D2)*4;

            min_of_3(D000,D2_cube[1][1][1],  Dpp0,D2_cube[2][2][1],  Dpm0,D2_cube[2][0][1]);
            du += (V_cube[2][2][1] - V_cube[2][0][1])*(D1-D2)*fx*fy*2;
            min_of_3(D000,D2_cube[1][1][1],  Dmp0,D2_cube[0][2][1],  Dmm0,D2_cube[0][0][1]);
            du += (V_cube[0][0][1] - V_cube[0][2][1])*(D1-D2)*fx*fy*2;
            //min_of_3(D000,D2_cube[1][1][1],  Dp0p,D2_cube[2][1][2],  Dp0m,D2_cube[2][1][0]);
            //du += (V_cube[2][1][2] - V_cube[2][1][0])*(D1-D2)*fx*fz*2;
            //min_of_3(D000,D2_cube[1][1][1],  Dm0p,D2_cube[0][1][2],  Dm0m,D2_cube[0][1][0]);
            //du += (V_cube[0][1][0] - V_cube[0][1][2])*(D1-D2)*fx*fz*2;
            //min_of_3(D000,D2_cube[1][1][1],  D0pp,D2_cube[1][2][2],  D0pm,D2_cube[1][2][0]);
            //du += (V_cube[1][2][2] - V_cube[1][2][0])*(D1-D2)*fy*fz*2;
            //min_of_3(D000,D2_cube[1][1][1],  D0mp,D2_cube[1][0][2],  D0mm,D2_cube[1][0][0]);
            //du += (V_cube[1][0][0] - V_cube[1][0][2])*(D1-D2)*fy*fz*2;
            min_of_3(D000,D2_cube[1][1][1],  Dp00,D2_cube[2][1][1],  Dm00,D2_cube[0][1][1]);
            du += (V_cube[2][1][1] - V_cube[0][1][1])*(D1-D2)*(2*fx*dfxdx+fy*dfxdy+fx*dfydy+ fz*dfxdz+fx*dfzdz);
            min_of_3(D000,D2_cube[1][1][1],  D0p0,D2_cube[1][2][1],  D0m0,D2_cube[1][0][1]);
            du += (V_cube[1][2][1] - V_cube[1][0][1])*(D1-D2)*(fx*dfydx+fy*dfxdx+2*fy*dfydy+ fz*dfydz+fy*dfzdz);
            //min_of_3(D000,D2_cube[1][1][1],  D00p,D2_cube[1][1][2],  D00m,D2_cube[1][1][0]);
            //du += (V_cube[1][1][2] - V_cube[1][1][0])*(D1-D2)*(fx*dfzdx+fz*dfxdx+fy*dfzdy+fz*dfydy+2*fz*dfzdz);
            //*/

            
            Courtemanche_Calculations(sparse_index);
            

            V_new_sparse[sparse_index] = V_cube[1][1][1] + (-I_ion/C_m + du/(4*dx*dx))*dt; // Transmembrane potential (mV)
            
            if((V_new_sparse[sparse_index] > 130.0) || (V_new_sparse[sparse_index] < -180.0) || (isinf(V_new_sparse[sparse_index])) || (isnan(V_new_sparse[sparse_index]))){
                printf("\nvoltage out of bound at x = %d, y = %d, z = %d, (V = %lf)\n", geo_x_index[sparse_index], geo_y_index[sparse_index], geo_z_index[sparse_index], V_new_sparse[sparse_index]);
                exit(1);
            }

        }
        }
        
        #pragma omp barrier
        for (sparse_index = 0; sparse_index < number_of_cells; sparse_index++){
            V_sparse[sparse_index] = V_new_sparse[sparse_index];
		}
        t += dt;
        increment++;
    }
    
    str = malloc(1000*sizeof(char));
    sprintf(str,"chmod -R 755 %s", main_output_directory_name);
    printf("%s\n",str);
    int statuss;
    statuss = system(str);
    if(save_to_mat_file == 1){
        char experiment_name[300];
        sprintf(experiment_name, "%s", c_file_name);
        char directory_to_read_from[300]; // the raw files that will be read
        char directory_to_save_mat_file_to[300]; // location to save the .mat file
        
        if(initial_conditions_from_files){
            sprintf(directory_to_read_from, "%s/L%d_continued_from_beat_%03d", main_run_directory, pacing_location, stimulus_to_continue);
        }
        else{
            sprintf(directory_to_read_from, "%s/L%d/",main_run_directory, pacing_location); //output
        }

        //sprintf(directory_to_save_mat_file_to, "%s/",main_run_directory);
        sprintf(directory_to_save_mat_file_to, "./matlab_variable_output");
        if (stat(directory_to_save_mat_file_to, &st) == -1) {
            mkdir(directory_to_save_mat_file_to, 0777);
        }

        char saved_mat_file_name[600];
        if(initial_conditions_from_files){
            sprintf(saved_mat_file_name, "%s_L%d_continued_from_beat_%03d.mat", experiment_name,pacing_location, stimulus_to_continue);
        }
        else{
            sprintf(saved_mat_file_name, "%s_L%d.mat", experiment_name, pacing_location);
        }
        
        char rel_path_to_mat_file[1000];
        sprintf(rel_path_to_mat_file, "%s/%s", directory_to_save_mat_file_to, saved_mat_file_name);
        
        str = malloc(1000*sizeof(char));
        sprintf(str, "matlab -batch \"non_void_reader_V3_auto_saver_function('%s','%s','%s','%s',%d)\"", geometry_file_name, directory_to_read_from, directory_to_save_mat_file_to, saved_mat_file_name, (int)round(dt*increment_divisor));
        printf("%s\n",str);
        printf("\n");
        int mat_file_run_tries = 0;
        int mat_file_exists = 0;
        while(1 == 1){
            statuss = system(str);
            in = fopen(rel_path_to_mat_file, "r");
            if(in != NULL){
                mat_file_exists = 1;
            }
            fclose(in);

            if((mat_file_exists == 1) || (mat_file_run_tries > 3)){
                if(mat_file_exists == 1){
                    sprintf(str, "chmod 755 %s",  rel_path_to_mat_file);
                    printf("%s\n", str);
                    statuss = system(str);
                }
                break;
            }
            else{
                mat_file_run_tries++;
                sleep(60);
            } 
        }
        if((delete_raw_files_after_saving_mat_file == 1) && (mat_file_exists == 1)){
            str = malloc(1000*sizeof(char));
            sprintf(str, "rm -r %s", directory_to_read_from);
            statuss = system(str);
        }
        
        if(make_figure_images_from_mat_file == 1){
            /*
            char directory_for_figure_images_folder[300]="saved_figure_images";
            if (stat(directory_for_figure_images_folder, &st) == -1) {
                mkdir(directory_for_figure_images_folder, 0777);
            }

            char images_folder_name[300];
            if(initial_conditions_from_files){
                sprintf(images_folder_name, "%s_L%d_continued_from_beat_%03d", experiment_name, pacing_location, stimulus_to_continue);
            }
            else{
                sprintf(images_folder_name, "%s_L%d", experiment_name, pacing_location);
            }

            char directory_for_figure_images[300];
            sprintf(directory_for_figure_images, "%s/%s", directory_for_figure_images_folder, images_folder_name);
            if (stat(directory_for_figure_images, &st) == -1) {
                mkdir(directory_for_figure_images, 0777);
            }

            sprintf(str, "matlab -batch \"save_figure_images('%s','%s','%s')\"", directory_to_save_mat_file_to, images_folder_name, directory_for_figure_images);
            printf("%s\n", str);
            statuss = system(str);
            */
        }
        free(str);
    }

    free(main_output_directory_name);
    free(c_file_name);
    free(saved_states_original_directory);
    free(saved_states_rerun_directory);
    free(main_run_directory);
    end_timer = time(NULL);
    printf("\ntime taken: %ld s\n", (end_timer - begin_timer));
    return 0;
}

void Courtemanche_Calculations (int sparse_index)
{   
    E_K  = ((R*T)/(F*z_K ))*log(     K_o/K_i_sparse[sparse_index]); // Equilibrium potential for K  (mV) // 131.867138
    E_Na = ((R*T)/(F*z_Na))*log(   Na_o/Na_i_sparse[sparse_index]); // Equilibrium potential for Na (mV)
    if(SAN_cell == 1){
        E_Ca = ((R*T)/(F*z_Ca))*log(Ca_o/Ca_i_SS_sparse[sparse_index]); // Equilibrium potential for Ca (mV)
    }
    else{
        E_Ca = ((R*T)/(F*z_Ca))*log(Ca_o/Ca_i_NJ_sparse[sparse_index]); // Equilibrium potential for Ca (mV)
    }
    
    Calculate_I_Na (sparse_index);
    //*
    Calculate_I_Kr (sparse_index);
    Calculate_I_Ks (sparse_index);
    Calculate_I_CaL (sparse_index);
    Calculate_I_CaT (sparse_index);
    Calculate_I_f (sparse_index);
    Calculate_I_to (sparse_index);
    Calculate_I_Kur (sparse_index);
    Calculate_I_K1 (sparse_index);
    Calculate_I_bCa (sparse_index);
    Calculate_I_bNa (sparse_index);
    Calculate_I_pCa (sparse_index);
    Calculate_I_NaCa (sparse_index);
    Calculate_I_NaK (sparse_index);
    Calculate_I_KACh (sparse_index);
    Calculate_Na_i_Concentration (sparse_index);
    Calculate_K_i_Concentration (sparse_index);
    Calculate_Ca_System_Concentrations (sparse_index);//*/
    Calculate_I_ion (sparse_index);
}

void Calculate_I_Na (int sparse_index) // Fast inward Na current (pA)
{
    if(SAN_cell == 1){
        // remember to covert from ms to s
        // remember to set the threadprivate variables
        // remember ICs
        // remember g_Na units and account for C_m being different for SAN (in here and final voltage)
        
        I_Na = scale_I_Na*C_m*g_Na_SAN*pow_3(m_sparse[sparse_index])*h_sparse[sparse_index]*(V_sparse[sparse_index] - E_mh_SAN);

        h_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+69.804)/4.4565));
        alpha_h = 20.0*exp(-0.125*(V_sparse[sparse_index]+75.0));
        beta_h = 2000.0/(320.0*exp(-0.1*(V_sparse[sparse_index]+75.0))+1.0);
        tau_h = 1.0/(alpha_h+beta_h);
        h_sparse[sparse_index] += (dt/1000.0)*(h_infinity - h_sparse[sparse_index])/tau_h;
        m_infinity = 1.0/(1.0+exp(-(V_sparse[sparse_index]+42.0504)/8.3106));
        if (fabs(V_sparse[sparse_index]+41.0) < 1.0e-5){alpha_m = 2000.0;}
        else{alpha_m = 200.0*(V_sparse[sparse_index]+41.0)/(1.0-exp(-0.1*(V_sparse[sparse_index]+41.0)));}
        beta_m = 8000.0*exp(-0.056*(V_sparse[sparse_index]+66.0));
        tau_m = 1.0/(alpha_m+beta_m);
        m_sparse[sparse_index] += (dt/1000.0)*(m_infinity-m_sparse[sparse_index])/tau_m;
    }
    else{

        I_Na = scale_I_Na*C_m*g_Na*pow_3(m_sparse[sparse_index])*h_sparse[sparse_index]*j_sparse[sparse_index]*(V_sparse[sparse_index] - E_Na);
        if(V_sparse[sparse_index] == -47.13){ alpha_m = 3.2;}
        else{ alpha_m = (0.32*(V_sparse[sparse_index] - INa_Act_Shift + 47.13)/(1 - exp(-0.1*(V_sparse[sparse_index] - INa_Act_Shift + 47.13))));}
        beta_m = 0.08*exp(-(V_sparse[sparse_index]-INa_Act_Shift)/11);
        m_infinity = alpha_m/(alpha_m+beta_m);
        tau_m = 1/(alpha_m+beta_m);
        m_sparse[sparse_index] = m_infinity - (m_infinity -  m_sparse[sparse_index])*exp(-dt/tau_m);
        if (V_sparse[sparse_index] < -40) {
            alpha_h = 0.135*exp(-(V_sparse[sparse_index] + 80.0)/6.8);
            beta_h = 3.56*exp(0.079*V_sparse[sparse_index]) + 310000*exp(0.35*V_sparse[sparse_index]);
        }
        else {
            alpha_h = 0.0;
            beta_h = 1.0 / (0.13*(1 + exp(-(V_sparse[sparse_index] + 10.66)/11.1)));
        }
        h_infinity = alpha_h/(alpha_h+beta_h);
        tau_h = 1/(alpha_h+beta_h);
        h_sparse[sparse_index] = h_infinity - (h_infinity - h_sparse[sparse_index])*exp(-dt/tau_h);
        if (V_sparse[sparse_index] < -40) {
            alpha_j = (-127140*exp(0.2444*V_sparse[sparse_index]) - 3.474e-5*exp(-0.04391*V_sparse[sparse_index])) * (V_sparse[sparse_index] + 37.78)/(1+exp(0.311*(V_sparse[sparse_index] + 79.23)));
            beta_j = 0.1212*exp(-0.01052*V_sparse[sparse_index]) / (1 + exp(-0.1378*(V_sparse[sparse_index] + 40.14)));
        }
        else {
            alpha_j = 0.0;
            beta_j = 0.3*exp(-2.535e-7*V_sparse[sparse_index]) / (1 + exp(-0.1*(V_sparse[sparse_index] + 32)));
        }
        j_infinity = alpha_j/(alpha_j+beta_j);
        tau_j = 1/(alpha_j+beta_j);
        j_sparse[sparse_index] = j_infinity - (j_infinity - j_sparse[sparse_index])*exp(-dt/tau_j);
    }
}

void Calculate_I_K1 (int sparse_index) // Inward rectifier K current (pA)
{
    if(SAN_cell == 1){
        I_K1 = scale_I_K1*C_m*g_K1_SAN*(V_sparse[sparse_index] - E_K_SAN) / (1 + exp(0.07*(V_sparse[sparse_index] + 80)));
    }
    else{
        I_K1 = scale_I_K1*C_m*g_K1*(V_sparse[sparse_index] - E_K + IK1_v_shift) / (1 + exp(0.07*(V_sparse[sparse_index] + 80 + IK1_v_shift)));
    }
}

void Calculate_I_to (int sparse_index) // Transient outward K current (pA)
{
    if(SAN_cell == 1){
        // remember to covert from ms to s
        // remember to set the threadprivate variables
        // remember ICs
        // remember g_Na units and account for C_m being different for SAN (in here and final voltage)
        // remember to use new g_Na_SAN variable
        
        I_to = scale_I_to*C_m*g_to_SAN*(V_sparse[sparse_index] - E_K_SAN)*oa_sparse[sparse_index]*oi_sparse[sparse_index];

        oa_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+49.0)/13.0));
        tau_oa = 0.001*0.6*(65.17/(0.57*exp(-0.08*(V_sparse[sparse_index]+44.0))+0.065*exp(0.1*(V_sparse[sparse_index]+45.93)))+10.1);   
        oa_sparse[sparse_index] += (dt/1000.0)*(oa_infinity-oa_sparse[sparse_index])/tau_oa;

        oi_infinity = 1.0/(1.0+exp(-(V_sparse[sparse_index]-19.3)/15.0));
        tau_oi = 0.001*0.66*1.4*(15.59/(1.037*exp(0.09*(V_sparse[sparse_index]+30.61))+0.369*exp(-0.12*(V_sparse[sparse_index]+23.84)))+2.98);
        oi_sparse[sparse_index] += (dt/1000.0)*(oi_infinity-oi_sparse[sparse_index])/tau_oi;
    }
    else{
        I_to = scale_I_to*1.05*C_m*g_to*oi_sparse[sparse_index]*oa_sparse[sparse_index]*(V_sparse[sparse_index] - E_K);
        
        oa_infinity = 1.0 / (1.0 + exp(((V_sparse[sparse_index] - Ito_Act_Shift) - 1.0) / -11.0));
        tau_oa = (0.0035 * exp(-((V_sparse[sparse_index] - Ito_Act_Shift) / 30.0) * 2) + 0.0015);
        oa_sparse[sparse_index] = oa_infinity + (oa_sparse[sparse_index] - oa_infinity) * exp(-(dt / 1000) / tau_oa);

        oi_infinity = 1.0 / (1.0 + exp((V_sparse[sparse_index] + 40.5) / 11.5));
        tau_oi = (0.025635 * exp(-((V_sparse[sparse_index] + 52.45) / 15.8827) * 2) + 0.01414);
        oi_sparse[sparse_index] = oi_infinity + (oi_sparse[sparse_index] - oi_infinity) * exp(-(dt / 1000) / tau_oi);
    }
}

void Calculate_I_Kur (int sparse_index)  // Ultrarapid delayed rectifier K current (pA)
{
    if(SAN_cell == 1){
        // remember to covert from ms to s
        // remember to set the threadprivate variables
        // remember ICs
        // remember g_Na units and account for C_m being different for SAN (in here and final voltage)
        // remember to use new g_Na_SAN variable
        
        I_Kur = scale_I_Kur*C_m*g_Kur_SAN*ua_sparse[sparse_index]*ui_sparse[sparse_index]*(V_sparse[sparse_index] - E_K_SAN); // 26=ua, 27=ui

        ua_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+6.0)/-8.6));
        tau_ua = 0.009/(1.0+exp((V_sparse[sparse_index]+5.0)/12.0))+0.0005;
        ua_sparse[sparse_index] += (dt/1000.0)*(ua_infinity-ua_sparse[sparse_index])/tau_ua;

        ui_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+7.5)/10.0));
        tau_ui = 0.59/(1.0+exp((V_sparse[sparse_index]+60.0)/10.0))+3.05;
		ui_sparse[sparse_index] += (dt/1000.0)*(ui_infinity-ui_sparse[sparse_index])/tau_ui;
    }
    else{
        I_Kur = scale_I_Kur*C_m*g_Kur*(4.5128 + 1.899769/(1.0+exp((V_sparse[sparse_index]-20.5232)/(-8.26597))))*ua_sparse[sparse_index]*ui_sparse[sparse_index]*(V_sparse[sparse_index] - E_K);

        ua_infinity = 1 / (1 + exp(-(V_sparse[sparse_index] + 5.52) / (8.6))); // change to simple version, haibo.  removed IKur mutations here.
        tau_ua = ((45.6666746826 / (1 + exp((V_sparse[sparse_index] + 11.2306497073) / 11.5254705962)) + 4.26753514993) * (0.262186042981 / (1 + exp((V_sparse[sparse_index] + 35.8658312707) / (-3.87510627762))) + 0.291755017928))/K_Q10;
        ua_sparse[sparse_index] = ua_infinity + (ua_sparse[sparse_index] - ua_infinity) * exp(-(dt) / tau_ua);
        
        ui_infinity = (0.52424) / (1.0 + exp((V_sparse[sparse_index] + 15.1142) / (7.567021))) + 0.4580778;
        tau_ui = (2328 / (1 + exp(((V_sparse[sparse_index]) - 9.435) / (3.5827))) + 1739.139)/K_Q10;
        ui_sparse[sparse_index] = ui_infinity + (ui_sparse[sparse_index] - ui_infinity) * exp(-(dt) / tau_ui);
    }
}


void Calculate_I_Kr (int sparse_index) // Rapid delayed rectifier K current (pA)
{
    if(SAN_cell == 1){
        I_Kr = scale_I_Kr*C_m*g_Kr_SAN*(V_sparse[sparse_index]-E_K_SAN)*(0.9*paF_sparse[sparse_index] + 0.1*paS_sparse[sparse_index])*xr_sparse[sparse_index];
        
        pa_infinity = 1.0/(1.0+exp(-(V_sparse[sparse_index]+10.0144)/7.6607));
        tau_paF = 1.0/(30.0*exp(V_sparse[sparse_index]/10.0)+exp(-V_sparse[sparse_index]/12.0));
        paF_sparse[sparse_index] += (dt/1000.0)*(pa_infinity-paF_sparse[sparse_index])/tau_paF; // paF (dimensionless) (in i_Kr_pa_gate)

        tau_paS = 0.84655354/(4.2*exp(V_sparse[sparse_index]/17.0)+0.15*exp(-V_sparse[sparse_index]/21.6));
        paS_sparse[sparse_index] += (dt/1000.0)*(pa_infinity-paS_sparse[sparse_index])/tau_paS; // paS (dimensionless) (in i_Kr_pa_gate)
        
        tau_xr = 1.0/(100.0*exp(-V_sparse[sparse_index]/54.645)+656.0*exp(V_sparse[sparse_index]/106.157));
        xr_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+28.6)/17.1));
        xr_sparse[sparse_index] += (dt/1000.0)*(xr_infinity-xr_sparse[sparse_index])/tau_xr; // piy (dimensionless) (in i_Kr_pi_gate)
    }
    else{
        I_Kr = scale_I_Kr*C_m *g_Kr*xr_sparse[sparse_index]*(V_sparse[sparse_index] - E_K) / (1+exp((V_sparse[sparse_index] + 15)/22.4));
        if(fabs(V_sparse[sparse_index] + 14.1) < 1e-10){alpha_xr = 0.0015;}
        else{alpha_xr = 0.0003*(V_sparse[sparse_index] + 14.1) / (1 - (exp(-(V_sparse[sparse_index] + 14.1)/5)));}
        if(fabs(V_sparse[sparse_index] - 3.3328) < 1e-10){beta_xr = 3.7836118e-4;}
        else{beta_xr = 7.3898e-5*(V_sparse[sparse_index]-3.3328)/(exp((V_sparse[sparse_index]-3.3328)/5.1237) - 1);}
        tau_xr = 1/(alpha_xr + beta_xr);
        xr_infinity = 1/(1 + exp(-(V_sparse[sparse_index] + 14.1)/6.5));
        xr_sparse[sparse_index] = xr_infinity-(xr_infinity-xr_sparse[sparse_index])*exp(-dt/tau_xr);
    }
}

void Calculate_I_Ks (int sparse_index) // Slow delayed rectifier K current (pA)
{
    if(SAN_cell == 1){
        I_Ks = scale_I_Ks*C_m*g_Ks_SAN*(V_sparse[sparse_index] - E_Ks_SAN)*pow_2(xs_sparse[sparse_index]); // 25=xs
        xs_infinity = sqrt(1.0/(1.0+exp(-(V_sparse[sparse_index]+0.6383)/10.7071)));
        alpha_xs = 28.0/(1.0+exp(-(V_sparse[sparse_index]-40.0)/3.0));
        beta_xs = 1.0*exp(-(V_sparse[sparse_index]-5.0)/25.0);
        tau_xs = 1.0/(alpha_xs+beta_xs);
		xs_sparse[sparse_index] += (dt/1000.0)*(xs_infinity-xs_sparse[sparse_index])/tau_xs;
    }
    else{
        I_Ks = scale_I_Ks*C_m*g_Ks*pow_2(xs_sparse[sparse_index])*(V_sparse[sparse_index] - E_K);
        if(fabs(V_sparse[sparse_index] - 19.9) < 1e-10){
            alpha_xs = 0.00068;
            beta_xs = 0.000315;
        }
        else{
            alpha_xs = 4e-5*(V_sparse[sparse_index]-19.9)/(1-exp(-(V_sparse[sparse_index]-19.9)/17));
            beta_xs = 3.5e-5*(V_sparse[sparse_index]-19.9)/(exp((V_sparse[sparse_index]-19.9)/9)-1);
        }
        tau_xs = 1/(2*(alpha_xs+beta_xs));
        xs_infinity = 1.0/sqrt(1 + exp(-(V_sparse[sparse_index]-19.9)/12.7));
        //xs_infinity = pow(1 + exp(-(V_sparse[sparse_index]-19.9)/12.7), -0.5);
        xs_sparse[sparse_index] = xs_infinity-(xs_infinity-xs_sparse[sparse_index])*exp(-dt/tau_xs);
    }
}

void Calculate_I_CaL (int sparse_index) // L-type inward Ca current (pA)
{
    if(SAN_cell == 1){
        //printf("%lf %lf %lf %lf %lf %lf\n", d_sparse[sparse_index],f_sparse[sparse_index],fCa_sparse[sparse_index],V_sparse[sparse_index],Ca_i_SS_sparse[sparse_index],exp(-2*V_sparse[sparse_index]*F/(R*T)));
        I_siCa = scale_I_CaL*C_m*(d_sparse[sparse_index]*f_sparse[sparse_index]*fCa_sparse[sparse_index])*2*(P_CaL*V_sparse[sparse_index]*F/(R*T))*(Ca_i_SS_sparse[sparse_index] - Ca_o*exp(-2*V_sparse[sparse_index]*F/(R*T)))/(1 - exp(-2*V_sparse[sparse_index]*F/(R*T)));
        I_siK = scale_I_CaL*C_m*(d_sparse[sparse_index]*f_sparse[sparse_index]*fCa_sparse[sparse_index])*(P_CaL*V_sparse[sparse_index]*F/(R*T))*0.000365*(K_i_SAN - K_o*exp(-V_sparse[sparse_index]*F/(R*T)))/(1 - exp(-V_sparse[sparse_index]*F/(R*T)));
        I_siNa = scale_I_CaL*C_m*(d_sparse[sparse_index]*f_sparse[sparse_index]*fCa_sparse[sparse_index])*(P_CaL*V_sparse[sparse_index]*F/(R*T))*0.0000185*(Na_i_SAN - Na_o*exp(-V_sparse[sparse_index]*F/(R*T)))/(1 - exp(-V_sparse[sparse_index]*F/(R*T)));
        I_CaL = (1.0 - ACh_block_I_CaL)*(I_siK+I_siNa+I_siCa);
        
        // d
        if (V_sparse[sparse_index] == -41.8){adVm = -41.80001;}
        else if (V_sparse[sparse_index] == -6.8){adVm = -6.80001;}
        else{adVm = V_sparse[sparse_index];}
        alpha_d = -0.02839*(adVm+41.8)/(exp(-(adVm+41.8)/2.5)-1.0)-0.0849*(adVm+6.8)/(exp(-(adVm+6.8)/4.8)-1.0);

        if (V_sparse[sparse_index] == -1.8){bdVm = -1.80001;}
        else{bdVm = V_sparse[sparse_index];}
        beta_d = 0.01143*(bdVm+1.8)/(exp((bdVm+1.8)/2.5)-1.0);

        d_infinity = 1.0/(1.0+exp(-(V_sparse[sparse_index]+16.4508)/4.3371));
        tau_d = 0.001/(alpha_d+beta_d);
        d_sparse[sparse_index] += (dt/1000.0)*(d_infinity-d_sparse[sparse_index])/tau_d; // dL (dimensionless) (in i_CaL_dL_gate)

        // f
        f_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+37.4)/5.3));
        tau_f = 0.001*(44.3+230.0*exp(-pow_2((V_sparse[sparse_index]+36.0)/10.0)));
        f_sparse[sparse_index] += (dt/1000.0)*(f_infinity-f_sparse[sparse_index])/tau_f;

        // fCa
        fCa_infinity = 0.000338/(0.000338 + Ca_i_SS_sparse[sparse_index]);
        tau_fCa = 0.001*fCa_infinity/0.0075;
        fCa_sparse[sparse_index] += (dt/1000.0)*(fCa_infinity-fCa_sparse[sparse_index])/tau_fCa;
    }
    else{
        I_CaL = scale_I_CaL*1.333333*C_m*g_CaL*d_sparse[sparse_index]*f_sparse[sparse_index]*fCa_sparse[sparse_index]*(V_sparse[sparse_index] - 65); //2017
        if(fabs(V_sparse[sparse_index] + 10) < 1e-10){tau_d = (4.579/(1 + exp(-(V_sparse[sparse_index] + 10)/6.24)));}
        else{tau_d = (1-exp(-(V_sparse[sparse_index] + 10)/6.24)) / (0.035*(V_sparse[sparse_index] + 10)*(1+exp(-(V_sparse[sparse_index] + 10)/6.24)));}
        d_infinity = 1/(1 + exp(-(V_sparse[sparse_index] + 10)/8));
        d_sparse[sparse_index] = d_infinity - (d_infinity - d_sparse[sparse_index])*exp(-dt/tau_d);

        tau_f = 9.0/(0.0197*exp(-pow_2(0.0337*(V_sparse[sparse_index] + 10))) + 0.02);
        f_infinity = 1/(1+exp((V_sparse[sparse_index] + 28)/6.9));
        f_sparse[sparse_index] = f_infinity - (f_infinity - f_sparse[sparse_index])*exp(-dt/tau_f);

        tau_fCa = 2.0;
        fCa_infinity = 1/(1 + Ca_i_SS_sparse[sparse_index]/0.00035); // 2017
        fCa_sparse[sparse_index] = fCa_infinity - (fCa_infinity - fCa_sparse[sparse_index])*exp(-dt/tau_fCa);
    }
}

void Calculate_I_CaT (int sparse_index){
    if(SAN_cell == 1){
        I_CaT = scale_I_CaT*2.0*C_m*P_CaT*(Ca_i_SS_sparse[sparse_index] - Ca_o*exp(-2.0*V_sparse[sparse_index]*F/(R*T))) *dT_sparse[sparse_index]*fT_sparse[sparse_index]* V_sparse[sparse_index]*F/(R*T*(1.0-exp(-1.0*V_sparse[sparse_index]*2.0*F/(R*T))));
        //I_CaT = 2.0*C_m*P_CaT*(1.0-Ca_o*exp(-2.0*V_sparse[sparse_index]*F/(R*T)))*dT_sparse[sparse_index]*fT_sparse[sparse_index]*V_sparse[sparse_index]*F/(R*T*(1.0-exp(-1.0*V_sparse[sparse_index]*2.0*F/(R*T))));
        
        dT_infinity = 1.0/(1.0+exp(-(V_sparse[sparse_index]+38.3)/5.5));
        tau_dT = 0.001/(1.068*exp((V_sparse[sparse_index]+38.3)/30.0)+1.068*exp(-(V_sparse[sparse_index]+38.3)/30.0));
        dT_sparse[sparse_index] += (dt/1000.0)*(dT_infinity-dT_sparse[sparse_index])/tau_dT; // dT (dimensionless) (in i_CaT_dT_gate)
        
        fT_infinity = 1.0/(1.0+exp((V_sparse[sparse_index]+58.7)/3.8));
        tau_fT = 1.0/(16.67*exp(-(V_sparse[sparse_index]+75.0)/83.3)+16.67*exp((V_sparse[sparse_index]+75.0)/15.38));
        fT_sparse[sparse_index] += (dt/1000.0)*(fT_infinity-fT_sparse[sparse_index])/tau_fT; // fT (dimensionless) (in i_CaT_fT_gate)
    }
    else{
        I_CaT = 0.0;
    }
}

void Calculate_I_f (int sparse_index){
    if(SAN_cell == 1){
        I_f_Na = scale_I_f*C_m*y_sparse[sparse_index]*g_f_Na*(V_sparse[sparse_index]-E_Na_SAN);
        I_f_K = scale_I_f*C_m*y_sparse[sparse_index]*g_f_K*(V_sparse[sparse_index]-E_K_SAN);
        I_f =  I_f_Na + I_f_K;
        tau_y = 1.0/(0.36*(V_sparse[sparse_index]+148.8-ACh_shift_I_f)/(exp(0.066*(V_sparse[sparse_index]+148.8-ACh_shift_I_f))-1.0)+0.1*(V_sparse[sparse_index]+87.3-ACh_shift_I_f)/(1.0-exp(-0.2*(V_sparse[sparse_index]+87.3-ACh_shift_I_f))))-0.054;
        if (V_sparse[sparse_index] < -(80.0-ACh_shift_I_f)){
            y_infinity = 0.01329+0.99921/(1.0+exp((V_sparse[sparse_index]+97.134-ACh_shift_I_f)/8.1752));
        }
        else{
            y_infinity = 0.0002501*exp(-(V_sparse[sparse_index]-ACh_shift_I_f)/12.861);
        }
        y_sparse[sparse_index] += (dt/1000.0)*(y_infinity-y_sparse[sparse_index])/tau_y;
    }
    else{
        I_f = 0.0;
    }
}

void Calculate_I_NaK (int sparse_index) // Na-K pump current (pA)
{
    if(SAN_cell == 1){
        //I_NaK = scale_I_NaK*C_m*I_NaKmax_SAN*pow(1.0+pow(K_mKo_SAN/K_o, 1.2), -1.0)*pow(1.0+pow(K_mNai_SAN/Na_i_SAN, 1.3), -1.0)*pow(1.0+exp(-(V_sparse[sparse_index]-E_Na_SAN+110.0)/20.0), -1.0);
        I_NaK = scale_I_NaK*C_m*I_NaKmax_SAN*(1.0/(1.0+pow(K_mKo_SAN/K_o, 1.2)))*(1.0/(1.0+pow(K_mNai_SAN/Na_i_SAN, 1.3)))*(1.0/(1.0+exp(-(V_sparse[sparse_index]-E_Na_SAN+110.0)/20.0)));
    }
    else{
        f_NaK = 1 / (1 + 0.1245*exp(-0.1*F*V_sparse[sparse_index]/(R*T)) + 0.0365*sigma*exp(-F*V_sparse[sparse_index]/(R*T)));
        I_NaK = scale_I_NaK*1.28*C_m*I_NaKmax*f_NaK*1/(1+pow_1p5(K_mNai/Na_i_sparse[sparse_index])) * K_o/(K_o+K_mKo); //2017
    }
}

void Calculate_I_NaCa (int sparse_index) // Na/Ca exchanger current (pA)
{
    if(SAN_cell == 1){
        k34 = Na_o/(4.663+Na_o);
        k32 = exp(0.4315*V_sparse[sparse_index]*F/(2.0*R*T));
        k43 = Na_i_SAN/(26.44+Na_i_SAN);
        di = 1.0+Ca_i_SS_sparse[sparse_index]/0.0207*(1.0+exp(-0.1369*V_sparse[sparse_index]*F/(R*T))+Na_i_SAN/26.44)+Na_i_SAN/395.3*(1.0+Na_i_SAN/2.289*(1.0+Na_i_SAN/26.44));
        k14 = Na_i_SAN/395.3*Na_i_SAN/2.289*(1.0+Na_i_SAN/26.44)*exp(0.4315*V_sparse[sparse_index]*F/(2.0*R*T))/di;
        k12 = Ca_i_SS_sparse[sparse_index]/0.0207*exp(-0.1369*V_sparse[sparse_index]*F/(R*T))/di;
        k41 = exp(-0.4315*V_sparse[sparse_index]*F/(2.0*R*T));
        x2 = k32*k43*(k14+k12) + k41*k12*(k34+k32);
        dodo = 1.0+Ca_o/3.663*2.0 + Na_o/1628.0*(1.0+Na_o/561.4*(1.0+Na_o/4.663));
        k21 = Ca_o/3.663*1.0/dodo;
        k23 = Na_o/1628.0*Na_o/561.4*(1.0+Na_o/4.663)*exp(-0.4315*V_sparse[sparse_index]*F/(2.0*R*T))/dodo;
        x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
        x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
        x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);
        I_NaCa = scale_I_NaCa*C_m*I_NaCamax_SAN*(x2*k21 - x1*k12)/(x1+x2+x3+x4);
    }
    else{
        I_NaCa = scale_I_NaCa*1.4*C_m*I_NaCamax*(exp(Gamma*F*V_sparse[sparse_index]/(R*T))*pow_3(Na_i_sparse[sparse_index])*Ca_o - exp((Gamma-1)*F*V_sparse[sparse_index]/(R*T))*pow_3(Na_o)*Ca_i_SS_sparse[sparse_index])
                                /  ((pow_3(K_mNa) + pow_3(Na_o))*(K_mCa+Ca_o)*(1 + k_sat*exp((Gamma-1)*F*V_sparse[sparse_index]/(R*T)))); //2017
    }
}

void Calculate_I_bCa (int sparse_index) // Background Ca current (pA)
{
    if(SAN_cell == 1){
        I_bCa = 0.0;
    }
    else{
        I_bCa = scale_I_bCa*C_m*g_bCa*(V_sparse[sparse_index] - E_Ca);
    }
}

void Calculate_I_bNa (int sparse_index) // Background Na current (pA)
{
    if(SAN_cell == 1){
        I_bNa = 0.0;
    }
    else{
        //if((geo_x_index[sparse_index]==126)&&(geo_y_index[sparse_index]==85)&&(geo_z_index[sparse_index]==2)){
        //    printf("%lf, %lf, %lf, %lf, %lf\n", scale_I_bNa, C_m, g_bNa, V_sparse[sparse_index], E_Na);
        //}
        I_bNa = scale_I_bNa*1.7*C_m*g_bNa*(V_sparse[sparse_index] - E_Na); //2017
    }
}

void Calculate_I_pCa (int sparse_index) // Sarcoplasmic Ca pump current (pA)
{
    if(SAN_cell == 1){
        I_pCa = 0.0;
    }
    else{
        I_pCa = scale_I_pCa*1.26*C_m*I_pCamax*Ca_i_SS_sparse[sparse_index]/(0.0005 + Ca_i_SS_sparse[sparse_index]); //2017
    }
}

void Calculate_I_KACh (int sparse_index) // Acetylcholine-activated K current (pA)
{
    if(fabs(ACh) > 1e-9){// from Human Atrial Action Potential and Ca2+ Model: Sinus Rhythm and Chronic Atrial Fibrillation, the model uses uM, in here we use mM
        if(SAN_cell == 1){
            I_KACh = scale_I_KACh*C_m*g_KACh_SAN*(V_sparse[sparse_index]-E_K_SAN)*(1.0+exp((V_sparse[sparse_index]+20.0)/20.0))*a_sparse[sparse_index];
            alpha_a = (3.5988-0.025641)/(1.0+0.0000012155/pow(1.0*ACh, 1.6951))+0.025641;
            beta_a = 10.0*exp(0.0133*(V_sparse[sparse_index]+40.0));
            a_infinity = alpha_a/(alpha_a+beta_a);
            tau_a = 1.0/(alpha_a+beta_a);
            a_sparse[sparse_index] += (dt/1000.0)*(a_infinity-a_sparse[sparse_index])/tau_a; // a (dimensionless) (in i_KACh_a_gate)
            
        }
        else{
            I_KACh = scale_I_KACh*C_m*(1/(1+pow(0.03/(ACh*1000.0), 2.1))) * (0.08 + 0.04/(1+exp((V_sparse[sparse_index]+91)/12))) * (V_sparse[sparse_index]-E_K);
        }
    }
    else{
        I_KACh = 0.0;
    }
}

void Calculate_I_ion (int sparse_index) // Total Ionic current (pA)
{
    if(SAN_cell == 1){
        I_ion = 1*I_Na+1*I_to+1*I_CaL+1*I_NaCa+1*I_Kur+1*I_Ks+1*I_NaK+1*I_Kr+1*I_KACh+1*I_f+1*I_CaT+1*I_K1; // temp
    }
    else{
        I_ion = 1*I_Na+1*I_to+1*I_CaL+1*I_NaCa+1*I_Kur+1*I_Ks+1*I_NaK+1*I_Kr+1*I_KACh+1*I_K1+1*I_pCa+1*I_bCa+1*I_bNa;    
    }
    
}

void Calculate_Na_i_Concentration (int sparse_index) // Intracellular concentration of Na (mM)
{
    if(SAN_cell == 1){
        Na_i_sparse[sparse_index] = 0.0;
    }
    else{
        Na_i_sparse[sparse_index] += -dt*(3*I_NaK + 3*I_NaCa + I_bNa + I_Na)/(F*Volume_i);
    }
}

void Calculate_K_i_Concentration (int sparse_index) // Intracellular concentration of K (mM)
{
    if(SAN_cell == 1){
        K_i_sparse[sparse_index] = 0.0;
    }
    else{
        K_i_sparse[sparse_index] +=  dt*(2*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks)/(F*Volume_i);
    }
}

void Calculate_Ca_System_Concentrations (int sparse_index)
{
    if(SAN_cell == 1){
        // change pow to pow_2
        j_SRCarel = 148041085.1*RyR_o_NJ_sparse[sparse_index]*(Ca_sr_SS_sparse[sparse_index]-Ca_i_SS_sparse[sparse_index]);
        diff = Ca_sr_SS_sparse[sparse_index]-Ca_i_SS_sparse[sparse_index];
        //kCaSR = 15.0-(15.0-1.0)/(1.0+ (pow(0.45/Ca_sr_SS_sparse[sparse_index], 2.5)) );
        kCaSR = 15.0-(15.0-1.0)/(1.0+ (sqrt(0.45/Ca_sr_SS_sparse[sparse_index])*pow_2(0.45/Ca_sr_SS_sparse[sparse_index])) );
        koSRCa = 10000.0/kCaSR;
        kiSRCa = 500.0*kCaSR;
        //delta_R = 5.0*RyR_a_NJ_sparse[sparse_index]-kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_a_SS_sparse[sparse_index]-(koSRCa*pow(Ca_i_SS_sparse[sparse_index], 2.0)*RyR_a_SS_sparse[sparse_index]-660.0*RyR_o_NJ_sparse[sparse_index]);
        //delta_O = koSRCa*pow(Ca_i_SS_sparse[sparse_index], 2.0)*RyR_a_SS_sparse[sparse_index]-660.0*RyR_o_NJ_sparse[sparse_index]-(kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_o_NJ_sparse[sparse_index]-5.0*RyR_o_SS_sparse[sparse_index]);
        //delta_I = kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_o_NJ_sparse[sparse_index]-5.0*RyR_o_SS_sparse[sparse_index]-(660.0*RyR_o_SS_sparse[sparse_index]-koSRCa*pow(Ca_i_SS_sparse[sparse_index], 2.0)*RyR_a_NJ_sparse[sparse_index]);
        //delta_RI = 660.0*RyR_o_SS_sparse[sparse_index]-koSRCa*pow(Ca_i_SS_sparse[sparse_index], 2.0)*RyR_a_NJ_sparse[sparse_index]-(5.0*RyR_a_NJ_sparse[sparse_index]-kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_a_SS_sparse[sparse_index]);
        delta_R = 5.0*RyR_a_NJ_sparse[sparse_index]-kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_a_SS_sparse[sparse_index]-(koSRCa*pow_2(Ca_i_SS_sparse[sparse_index])*RyR_a_SS_sparse[sparse_index]-660.0*RyR_o_NJ_sparse[sparse_index]);
        delta_O = koSRCa*pow_2(Ca_i_SS_sparse[sparse_index])*RyR_a_SS_sparse[sparse_index]-660.0*RyR_o_NJ_sparse[sparse_index]-(kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_o_NJ_sparse[sparse_index]-5.0*RyR_o_SS_sparse[sparse_index]);
        delta_I = kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_o_NJ_sparse[sparse_index]-5.0*RyR_o_SS_sparse[sparse_index]-(660.0*RyR_o_SS_sparse[sparse_index]-koSRCa*pow_2(Ca_i_SS_sparse[sparse_index])*RyR_a_NJ_sparse[sparse_index]);
        delta_RI = 660.0*RyR_o_SS_sparse[sparse_index]-koSRCa*pow_2(Ca_i_SS_sparse[sparse_index])*RyR_a_NJ_sparse[sparse_index]-(5.0*RyR_a_NJ_sparse[sparse_index]-kiSRCa*Ca_i_SS_sparse[sparse_index]*RyR_a_SS_sparse[sparse_index]);
        P_tot = RyR_o_SS_sparse[sparse_index]+RyR_o_NJ_sparse[sparse_index]+RyR_a_SS_sparse[sparse_index]+RyR_a_NJ_sparse[sparse_index];
        delta_fTC = 88800.0*Ca_i_NJ_sparse[sparse_index]*(1.0-fTC_sparse[sparse_index])-446.0*fTC_sparse[sparse_index];
        delta_fTMC = 227700.0*Ca_i_NJ_sparse[sparse_index]*(1.0-(Ca_serca_SS_sparse[sparse_index]+Ca_serca_NJ_sparse[sparse_index]))-7.51*Ca_serca_SS_sparse[sparse_index];
        delta_fTMM = 2277.0*2.5*(1.0-(Ca_serca_SS_sparse[sparse_index]+Ca_serca_NJ_sparse[sparse_index]))-751.0*Ca_serca_NJ_sparse[sparse_index];
        delta_fCMi = 1.642e6*Ca_i_NJ_sparse[sparse_index]*(1.0-RyR_c_SS_sparse[sparse_index])-542.0*RyR_c_SS_sparse[sparse_index];

        delta_fCMs = 1.642e6*Ca_i_SS_sparse[sparse_index]*(1.0-RyR_c_NJ_sparse[sparse_index]) - 542.0*RyR_c_NJ_sparse[sparse_index];
        //printf("%lf = %lf*%lf - %lf\n", delta_fCMs, 1.642e6*Ca_i_NJ_sparse[sparse_index],1.0-RyR_c_SS_sparse[sparse_index], 542.0*RyR_c_SS_sparse[sparse_index]);

        delta_fCQ = 175.4*Ca_sr_SS_sparse[sparse_index]*(1.0-fCQ_sparse[sparse_index])-445.0*fCQ_sparse[sparse_index];
        j_Ca_dif = (Ca_i_SS_sparse[sparse_index]-Ca_i_NJ_sparse[sparse_index])/5.469e-5;
        if(fabs(ACh) > 1e-9){
            ACh_block_P_up = (0.7*ACh/(0.00009+ACh));
        }
        else{
            ACh_block_P_up = 0.0;
        }
        j_up = (5.0*(1.0- ACh_block_P_up))/(1.0+exp((-Ca_i_NJ_sparse[sparse_index]+0.000286113)/5.0e-5));
        delta_Ca_i = 1.0*(j_Ca_dif*V_sub-j_up*V_nsr)/V_i-(0.045*delta_fCMi+0.031*delta_fTC+0.062*delta_fTMC);
        
        delta_Ca_sub = j_SRCarel*V_jsr/V_sub-(1e-6*(I_siCa+I_CaT-2.0*I_NaCa)/(2.0*F*V_sub)+ j_Ca_dif+0.045*delta_fCMs);
        //delta_Ca_sub = 0.3747*V_jsr/V_sub-(1e-6*(-4.1658 + -1.5079 -2.0*-4.863)/(2.0*F*V_sub)+ 0.9710 +0.045*0.108);
        //printf("%lf = %lf %lf %lf %lf %lf %lf\n", delta_Ca_sub, j_SRCarel, I_siCa, I_CaT, I_NaCa, j_Ca_dif, delta_fCMs);

        j_tr = (Ca_sr_NJ_sparse[sparse_index]-Ca_sr_SS_sparse[sparse_index])/0.04;
        delta_Ca_nsr = j_up-j_tr*V_jsr/V_nsr;
        delta_Ca_jsr = j_tr-(j_SRCarel+10.0*delta_fCQ);

        //*
        RyR_o_SS_sparse[sparse_index]  +=  (dt/1000.0)*delta_I; // I (dimensionless) (in Ca_SR_release)
        RyR_o_NJ_sparse[sparse_index]  +=  (dt/1000.0)*delta_O; // O (dimensionless) (in Ca_SR_release)
        RyR_a_SS_sparse[sparse_index]  +=  (dt/1000.0)*delta_R; // R1 (dimensionless) (R in Ca_SR_release)
        RyR_a_NJ_sparse[sparse_index]  +=  (dt/1000.0)*delta_RI; // RI (dimensionless) (in Ca_SR_release)
        RyR_c_SS_sparse[sparse_index]  +=  (dt/1000.0)*delta_fCMi; // fCMi (dimensionless) (in Ca_buffering)
        RyR_c_NJ_sparse[sparse_index]  +=  (dt/1000.0)*delta_fCMs; // fCMs (dimensionless) (in Ca_buffering)
        fCQ_sparse[sparse_index]  +=  (dt/1000.0)*delta_fCQ; // fCQ (dimensionless) (in Ca_buffering)
        fTC_sparse[sparse_index]  +=  (dt/1000.0)*delta_fTC; // fTC (dimensionless) (in Ca_buffering)
        Ca_serca_SS_sparse[sparse_index]  +=  (dt/1000.0)*delta_fTMC; // fTMC (dimensionless) (in Ca_buffering)
        Ca_serca_NJ_sparse[sparse_index]  +=  (dt/1000.0)*delta_fTMM; // fTMM (dimensionless) (in Ca_buffering)
        Ca_sr_SS_sparse[sparse_index] +=  (dt/1000.0)*delta_Ca_jsr; // Ca_jsr (millimolar) (in Ca_dynamics)
        Ca_sr_NJ_sparse[sparse_index] +=  (dt/1000.0)*delta_Ca_nsr; // Ca_nsr (millimolar) (in Ca_dynamics)

        Ca_i_SS_sparse[sparse_index] +=  (dt/1000.0)*delta_Ca_sub; // Ca_sub (millimolar) (in Ca_dynamics), the only variable that becomes unstable at large dt

        Ca_i_NJ_sparse[sparse_index] +=  (dt/1000.0)*delta_Ca_i; // Ca_i (millimolar) (in Ca_dynamics)
        //*/
    }
    else{
        BSS = 1.0/( 1 + SLlow*KdSLlow/pow_2(Ca_i_SS_sparse[sparse_index]+KdSLlow) + SLhigh*KdSLhigh/pow_2(Ca_i_SS_sparse[sparse_index]+KdSLhigh) + BCa*KdBCa/pow_2(Ca_i_SS_sparse[sparse_index]+KdBCa));
        B_i_NJ = 1.0/(1 + BCa*KdBCa/pow_2(Ca_i_NJ_sparse[sparse_index]+KdBCa));
        B_sr_NJ = 1.0/(1 + CSQN*KdCSQN/pow_2(Ca_sr_NJ_sparse[sparse_index]+ KdCSQN));
        B_sr_SS = 1.0/(1 + CSQN*KdCSQN/pow_2(Ca_sr_SS_sparse[sparse_index] + KdCSQN));
        I_serca_SR_NJ = SERCA_Scale*1.5*(-k3*pow_2(Ca_sr_NJ_sparse[sparse_index])*(Cpumps-Ca_serca_NJ_sparse[sparse_index]) + k4*Ca_serca_NJ_sparse[sparse_index])*V_NJ;
        I_serca_bulk_NJ = SERCA_Scale*1.5*(k1*pow_2(Ca_i_NJ_sparse[sparse_index])*(Cpumps-Ca_serca_NJ_sparse[sparse_index]) - k2*Ca_serca_NJ_sparse[sparse_index])*V_NJ;
        I_serca_SR_SS = SERCA_Scale*1.5*(-k3*pow_2(Ca_sr_SS_sparse[sparse_index])*(Cpumps-Ca_serca_SS_sparse[sparse_index]) + k4*Ca_serca_SS_sparse[sparse_index])*V_SS;
        I_serca_bulk_SS = SERCA_Scale*1.5*(k1*pow_2(Ca_i_SS_sparse[sparse_index])*(Cpumps-Ca_serca_SS_sparse[sparse_index]) - k2*Ca_serca_SS_sparse[sparse_index])*V_SS;
        I_SS_NJ = 2.5*DCa*ASS_NJ/XSS_NJ*((Ca_i_SS_sparse[sparse_index])-Ca_i_NJ_sparse[sparse_index])*1e-6;
        RyR_serca_SS = 1 - 1/(1 + exp((Ca_sr_SS_sparse[sparse_index] - 0.3) / 0.1));
        RyR_a_SS_infinity = 0.505 - 0.427/( 1 + exp(((1.0 + fRyR)*Ca_i_SS_sparse[sparse_index]*1000 - 0.29)/0.082) );
        RyR_o_SS_infinity = 1 - 1/(1 + exp(((1.0 + fRyR)*Ca_i_SS_sparse[sparse_index]*1000 - RyR_a_SS_sparse[sparse_index] - 0.22)/0.03));
        RyR_c_SS_infinity = 1 / (1 + exp(((1.0 + fRyR)*Ca_i_SS_sparse[sparse_index]*1000 - RyR_a_SS_sparse[sparse_index] - 0.02)/0.01));  
        I_rel_SS = fIRel*625*V_SS*RyR_o_SS_sparse[sparse_index]*RyR_c_SS_sparse[sparse_index]*RyR_serca_SS*(Ca_sr_SS_sparse[sparse_index] - Ca_i_SS_sparse[sparse_index]);
        RyR_serca_NJ = 1 - 1/(1 + exp((Ca_sr_NJ_sparse[sparse_index] - 0.3)/0.1));
        RyR_a_NJ_infinity = 0.505 - 0.427/(1 + exp(((2 + fRyR)*Ca_i_NJ_sparse[sparse_index]*1000 - 0.29)/0.082));
        RyR_o_NJ_infinity = 1 - 1/(1 + exp(((2 + fRyR)*Ca_i_NJ_sparse[sparse_index]*1000 - RyR_a_NJ_sparse[sparse_index] - 0.22)/0.03));
        RyR_c_NJ_infinity = 1 / (1 + exp(((2 + fRyR)*Ca_i_NJ_sparse[sparse_index]*1000 - RyR_a_NJ_sparse[sparse_index] - 0.02)/0.01));
        I_rel_NJ = fIRel*V_NJ*RyR_o_NJ_sparse[sparse_index]*RyR_c_NJ_sparse[sparse_index]*RyR_serca_NJ*(Ca_sr_NJ_sparse[sparse_index] - Ca_i_NJ_sparse[sparse_index]);
        I_SRleak_NJ = GSR_leak*0.5*kSRleak * (Ca_sr_NJ_sparse[sparse_index] - Ca_i_NJ_sparse[sparse_index]) * V_NJ;
        I_SRleak_SS = GSR_leak*0.5*kSRleak * (Ca_sr_SS_sparse[sparse_index] - Ca_i_SS_sparse[sparse_index]) * V_SS;
        I_Ca_NJ = -I_serca_bulk_NJ + I_SRleak_NJ + I_rel_NJ + I_SS_NJ;
        I_Ca_SS = -I_SS_NJ + I_SRleak_SS - I_serca_bulk_SS + I_rel_SS;
        I_SRCa_NJ = I_serca_SR_NJ - I_SRleak_NJ - I_rel_NJ;
        I_SRCa_SS = I_serca_SR_SS - I_SRleak_SS - I_rel_SS;
        Ca_serca_SS_sparse[sparse_index] += (dt / 1000) * (I_serca_bulk_SS - I_serca_SR_SS)/V_SS;;
        Ca_serca_NJ_sparse[sparse_index] += (dt / 1000) * (I_serca_bulk_NJ - I_serca_SR_NJ)/V_NJ;
        RyR_a_SS_sparse[sparse_index] = RyR_a_SS_infinity - (RyR_a_SS_infinity - RyR_a_SS_sparse[sparse_index]) * exp(-(dt / 1000) / (TRyRadapt*RyRTauScale));
        RyR_a_NJ_sparse[sparse_index] = RyR_a_NJ_infinity - (RyR_a_NJ_infinity - RyR_a_NJ_sparse[sparse_index]) * exp(-(dt / 1000) / (TRyRadapt*RyRTauScale));
        RyR_o_SS_sparse[sparse_index] = RyR_o_SS_infinity - (RyR_o_SS_infinity - RyR_o_SS_sparse[sparse_index]) * exp(-(dt / 1000) / TRyRactSS);
        RyR_o_NJ_sparse[sparse_index] = RyR_o_NJ_infinity - (RyR_o_NJ_infinity - RyR_o_NJ_sparse[sparse_index]) * exp(-(dt / 1000) / TRyRactNJ);
        RyR_c_SS_sparse[sparse_index] = RyR_c_SS_infinity - (RyR_c_SS_infinity - RyR_c_SS_sparse[sparse_index]) * exp(-(dt / 1000) / (TRyRinactSS*RyRTauScale));
        RyR_c_NJ_sparse[sparse_index] = RyR_c_NJ_infinity - (RyR_c_NJ_infinity - RyR_c_NJ_sparse[sparse_index]) * exp(-(dt / 1000) / (TRyRinactNJ*RyRTauScale));
        Ca_i_SS_sparse[sparse_index] += (dt / 1000.0) * BSS * (I_Ca_SS/V_SS + (-I_CaL-I_bCa-I_pCa+2*I_NaCa) / (2*V_SS*1000*F));  ;
        Ca_i_NJ_sparse[sparse_index] += (dt / 1000.0) * I_Ca_NJ / V_NJ * B_i_NJ;
        dCa_SR_NJ = dt*( 5*B_sr_NJ*DCaSR*(Ca_sr_SS_sparse[sparse_index] - Ca_sr_NJ_sparse[sparse_index])/(6*r*r) + B_sr_NJ*I_SRCa_NJ/V_SR_NJ);
        dCa_SR_SS = dt*( 7*B_sr_SS*DCaSR*(Ca_sr_NJ_sparse[sparse_index] - Ca_sr_SS_sparse[sparse_index])/(8*r*r) + B_sr_SS*I_SRCa_SS/V_SR_SS);
        Ca_sr_NJ_sparse[sparse_index] += dCa_SR_NJ;
        Ca_sr_SS_sparse[sparse_index] += dCa_SR_SS;
    }
}
