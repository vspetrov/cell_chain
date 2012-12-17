//==============================================================================
// CellML file:   C:\Users\valentin\Downloads\maltsev_2009_steady_state.cellml
// CellML model:  maltsev_2009
// Date and time: 28.08.2012 at 14:14:18
//------------------------------------------------------------------------------
// Conversion from CellML 1.0 to C++ was done using COR (0.9.31.1409)
//    Copyright 2002-2012 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================

#include "maltsev_2009.h"

//------------------------------------------------------------------------------

#include <math.h>

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void maltsev_2009::init()
{
   //---------------------------------------------------------------------------
   // State variables
   //---------------------------------------------------------------------------

   Y[0] = 0.694241313965;   // q (dimensionless) (in AP_sensitive_currents_q_gate)
   Y[1] = 0.00558131733359;   // r (dimensionless) (in AP_sensitive_currents_r_gate)
   Y[2] = -57.9639346865;   // Vm (millivolt) (in Vm)
   Y[3] = 0.0594880901438;   // fCMi (dimensionless) (in calcium_buffering)
   Y[4] = 0.054381370046;   // fCMs (dimensionless) (in calcium_buffering)
   Y[5] = 0.273207128393;   // fCQ (dimensionless) (in calcium_buffering)
   Y[6] = 0.0291316176172;   // fTC (dimensionless) (in calcium_buffering)
   Y[7] = 0.432694959597;   // fTMC (dimensionless) (in calcium_buffering)
   Y[8] = 0.501049376634;   // fTMM (dimensionless) (in calcium_buffering)
   Y[9] = 0.316762674605;   // Ca_jsr (millimolar) (in calcium_dynamics)
   Y[10] = 1.49348117734;   // Ca_nsr (millimolar) (in calcium_dynamics)
   Y[11] = 0.000138112560112;   // Ca_sub (millimolar) (in calcium_dynamics)
   Y[12] = 0.000150018670943;   // Cai (millimolar) (in calcium_dynamics)
   Y[13] = 0.000584545564405;   // dL (dimensionless) (in i_CaL_dL_gate)
   Y[14] = 0.711395919653;   // fCa (dimensionless) (in i_CaL_fCa_gate)
   Y[15] = 0.862381249774;   // fL (dimensionless) (in i_CaL_fL_gate)
   Y[16] = 0.00504393374639;   // dT (dimensionless) (in i_CaT_dT_gate)
   Y[17] = 0.420757825415;   // fT (dimensionless) (in i_CaT_fT_gate)
   Y[18] = 0.144755091176;   // paF (dimensionless) (in i_Kr_pa_gate)
   Y[19] = 0.453100576739;   // paS (dimensionless) (in i_Kr_pa_gate)
   Y[20] = 0.849409822329;   // pi_ (dimensionless) (in i_Kr_pi_gate)
   Y[21] = 0.0264600410928;   // n (dimensionless) (in i_Ks_n_gate)
   Y[22] = 0.113643187247;   // y (dimensionless) (in i_f_y_gate)
   Y[23] = 0.42380243163;   // qa (dimensionless) (in i_st_qa_gate)
   Y[24] = 0.447294008304;   // qi (dimensionless) (in i_st_qi_gate)
   Y[25] = 7.86181717518e-8;   // I (dimensionless) (in j_SRCarel)
   Y[26] = 1.7340201253e-7;   // O (dimensionless) (in j_SRCarel)
   Y[27] = 0.688047760973;   // R1 (dimensionless) (R in j_SRCarel)
   Y[28] = 0.311951987007;   // RI (dimensionless) (in j_SRCarel)

   YNames[0].assign("q");
   YNames[1].assign("r");
   YNames[2].assign("Vm");
   YNames[3].assign("fCMi");
   YNames[4].assign("fCMs");
   YNames[5].assign("fCQ");
   YNames[6].assign("fTC");
   YNames[7].assign("fTMC");
   YNames[8].assign("fTMM");
   YNames[9].assign("Ca_jsr");
   YNames[10].assign("Ca_nsr");
   YNames[11].assign("Ca_sub");
   YNames[12].assign("Cai");
   YNames[13].assign("dL");
   YNames[14].assign("fCa");
   YNames[15].assign("fL");
   YNames[16].assign("dT");
   YNames[17].assign("fT");
   YNames[18].assign("paF");
   YNames[19].assign("paS");
   YNames[20].assign("pi_");
   YNames[21].assign("n");
   YNames[22].assign("y");
   YNames[23].assign("qa");
   YNames[24].assign("qi");
   YNames[25].assign("I");
   YNames[26].assign("O");
   YNames[27].assign("R1");
   YNames[28].assign("RI");

   YUnits[0].assign("dimensionless");
   YUnits[1].assign("dimensionless");
   YUnits[2].assign("millivolt");
   YUnits[3].assign("dimensionless");
   YUnits[4].assign("dimensionless");
   YUnits[5].assign("dimensionless");
   YUnits[6].assign("dimensionless");
   YUnits[7].assign("dimensionless");
   YUnits[8].assign("dimensionless");
   YUnits[9].assign("millimolar");
   YUnits[10].assign("millimolar");
   YUnits[11].assign("millimolar");
   YUnits[12].assign("millimolar");
   YUnits[13].assign("dimensionless");
   YUnits[14].assign("dimensionless");
   YUnits[15].assign("dimensionless");
   YUnits[16].assign("dimensionless");
   YUnits[17].assign("dimensionless");
   YUnits[18].assign("dimensionless");
   YUnits[19].assign("dimensionless");
   YUnits[20].assign("dimensionless");
   YUnits[21].assign("dimensionless");
   YUnits[22].assign("dimensionless");
   YUnits[23].assign("dimensionless");
   YUnits[24].assign("dimensionless");
   YUnits[25].assign("dimensionless");
   YUnits[26].assign("dimensionless");
   YUnits[27].assign("dimensionless");
   YUnits[28].assign("dimensionless");

   YComponents[0].assign("AP_sensitive_currents_q_gate");
   YComponents[1].assign("AP_sensitive_currents_r_gate");
   YComponents[2].assign("Vm");
   YComponents[3].assign("calcium_buffering");
   YComponents[4].assign("calcium_buffering");
   YComponents[5].assign("calcium_buffering");
   YComponents[6].assign("calcium_buffering");
   YComponents[7].assign("calcium_buffering");
   YComponents[8].assign("calcium_buffering");
   YComponents[9].assign("calcium_dynamics");
   YComponents[10].assign("calcium_dynamics");
   YComponents[11].assign("calcium_dynamics");
   YComponents[12].assign("calcium_dynamics");
   YComponents[13].assign("i_CaL_dL_gate");
   YComponents[14].assign("i_CaL_fCa_gate");
   YComponents[15].assign("i_CaL_fL_gate");
   YComponents[16].assign("i_CaT_dT_gate");
   YComponents[17].assign("i_CaT_fT_gate");
   YComponents[18].assign("i_Kr_pa_gate");
   YComponents[19].assign("i_Kr_pa_gate");
   YComponents[20].assign("i_Kr_pi_gate");
   YComponents[21].assign("i_Ks_n_gate");
   YComponents[22].assign("i_f_y_gate");
   YComponents[23].assign("i_st_qa_gate");
   YComponents[24].assign("i_st_qi_gate");
   YComponents[25].assign("j_SRCarel");
   YComponents[26].assign("j_SRCarel");
   YComponents[27].assign("j_SRCarel");
   YComponents[28].assign("j_SRCarel");

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   g_sus = 0.02;   // nanoS_per_picoF (in AP_sensitive_currents)
   g_to = 0.252;   // nanoS_per_picoF (in AP_sensitive_currents)
   CM_tot = 0.045;   // millimolar (in calcium_buffering)
   CQ_tot = 10.0;   // millimolar (in calcium_buffering)
   TC_tot = 0.031;   // millimolar (in calcium_buffering)
   TMC_tot = 0.062;   // millimolar (in calcium_buffering)
   kb_CM = 0.542;   // per_millisecond (in calcium_buffering)
   kb_CQ = 0.445;   // per_millisecond (in calcium_buffering)
   kb_TC = 0.446;   // per_millisecond (in calcium_buffering)
   kb_TMC = 0.00751;   // per_millisecond (in calcium_buffering)
   kb_TMM = 0.751;   // per_millisecond (in calcium_buffering)
   kf_CM = 227.7;   // per_millimolar_millisecond (in calcium_buffering)
   kf_CQ = 0.534;   // per_millimolar_millisecond (in calcium_buffering)
   kf_TC = 88.8;   // per_millimolar_millisecond (in calcium_buffering)
   kf_TMC = 227.7;   // per_millimolar_millisecond (in calcium_buffering)
   kf_TMM = 2.277;   // per_millimolar_millisecond (in calcium_buffering)
   Km_fCa = 0.00035;   // millimolar (in i_CaL_fCa_gate)
   alpha_fCa = 0.021;   // per_millisecond (in i_CaL_fCa_gate)
   E_CaL = 45.0;   // millivolt (in i_CaL)
   g_CaL = 0.464;   // nanoS_per_picoF (in i_CaL)
   E_CaT = 45.0;   // millivolt (in i_CaT)
   g_CaT = 0.1832;   // nanoS_per_picoF (in i_CaT)
   g_Kr = 0.08113973;   // nanoS_per_picoF (in i_Kr)
   g_Ks = 0.0259;   // nanoS_per_picoF (in i_Ks)
   K1ni = 395.3;   // millimolar (in i_NaCa)
   K1no = 1628.0;   // millimolar (in i_NaCa)
   K2ni = 2.289;   // millimolar (in i_NaCa)
   K2no = 561.4;   // millimolar (in i_NaCa)
   K3ni = 26.44;   // millimolar (in i_NaCa)
   K3no = 4.663;   // millimolar (in i_NaCa)
   Kci = 0.0207;   // millimolar (in i_NaCa)
   Kcni = 26.44;   // millimolar (in i_NaCa)
   Kco = 3.663;   // millimolar (in i_NaCa)
   Qci = 0.1369;   // dimensionless (in i_NaCa)
   Qco = 0.0;   // dimensionless (in i_NaCa)
   Qn = 0.4315;   // dimensionless (in i_NaCa)
   kNaCa = 187.5;   // picoA_per_picoF (in i_NaCa)
   Km_Kp = 1.4;   // millimolar (in i_NaK)
   Km_Nap = 14.0;   // millimolar (in i_NaK)
   i_NaK_max = 2.88;   // picoA_per_picoF (in i_NaK)
   g_b_Ca = 0.0006;   // nanoS_per_picoF (in i_b_Ca)
   g_b_Na = 0.00486;   // nanoS_per_picoF (in i_b_Na)
   VIf_half = -64.0;   // millivolt (in i_f_y_gate)
   g_if = 0.15;   // nanoS_per_picoF (in i_f)
   E_st = 37.4;   // millivolt (in i_st)
   g_st = 0.003;   // nanoS_per_picoF (in i_st)
   K_up = 0.0006;   // millimolar (in intracellular_calcium_fluxes)
   P_up = 0.012;   // millimolar_per_millisecond (in intracellular_calcium_fluxes)
   tau_dif_Ca = 0.04;   // millisecond (in intracellular_calcium_fluxes)
   tau_tr = 40.0;   // millisecond (in intracellular_calcium_fluxes)
   EC50_SR = 0.45;   // millimolar (in j_SRCarel)
   HSR = 2.5;   // dimensionless (in j_SRCarel)
   MaxSR = 15.0;   // dimensionless (in j_SRCarel)
   MinSR = 1.0;   // dimensionless (in j_SRCarel)
   kiCa = 0.5;   // per_millimolar_millisecond (in j_SRCarel)
   kim = 0.005;   // per_millisecond (in j_SRCarel)
   koCa = 10.0;   // per_millimolar2_millisecond (in j_SRCarel)
   kom = 0.06;   // per_millisecond (in j_SRCarel)
   ks = 250000.0;   // per_millisecond (in j_SRCarel)
   Cao = 2.0;   // millimolar (in model_parameters)
   Cm = 32.0;   // picoF (in model_parameters)
   F = 96485.0;   // coulomb_per_mole (in model_parameters)
   Ki = 140.0;   // millimolar (in model_parameters)
   Ko = 5.4;   // millimolar (in model_parameters)
   L_cell = 70.0;   // micrometre (in model_parameters)
   L_sub = 0.02;   // micrometre (in model_parameters)
   Mgi = 2.5;   // millimolar (in model_parameters)
   Nai = 10.0;   // millimolar (in model_parameters)
   Nao = 140.0;   // millimolar (in model_parameters)
   R2 = 8314.4;   // joule_per_kilomole_kelvin (R in model_parameters)
   R_cell = 4.0;   // micrometre (in model_parameters)
   T = 310.15;   // kelvin (in model_parameters)
   V_i_part = 0.46;   // dimensionless (in model_parameters)
   V_jsr_part = 0.0012;   // dimensionless (in model_parameters)
   V_nsr_part = 0.0116;   // dimensionless (in model_parameters)

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   E_K = R2*T/F*log(Ko/Ki);
   E_Na = R2*T/F*log(Nao/Nai);
   E_Ks = R2*T/F*log((Ko+0.12*Nao)/(Ki+0.12*Nai));
   RTOnF = R2*T/F;
   k43 = Nai/(K3ni+Nai);
   k34 = Nao/(K3no+Nao);
   V_sub = 0.001*2.0*3.14159265358979*L_sub*(R_cell-L_sub/2.0)*L_cell;
   V_cell = 0.001*3.14159265358979*pow(R_cell, 2.0)*L_cell;
   V_nsr = V_nsr_part*V_cell;
   V_i = V_i_part*V_cell-V_sub;
   V_jsr = V_jsr_part*V_cell;
}

//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

void maltsev_2009::compute(double time)
{
   // time: time (millisecond)

   i_to = Cm*g_to*(Y[2]-E_K)*Y[0]*Y[1];
   i_sus = Cm*g_sus*(Y[2]-E_K)*Y[1];
   q_infinity = 1.0/(1.0+exp((Y[2]+49.0)/13.0));
   tau_q = 6.06+39.102/(0.57*exp(-0.08*(Y[2]+44.0))+0.065*exp(0.1*(Y[2]+45.93)));
   dY[0] = (q_infinity-Y[0])/tau_q;
   r_infinity = 1.0/(1.0+exp(-(Y[2]-19.3)/15.0));
   tau_r = 2.75352+14.40516/(1.037*exp(0.09*(Y[2]+30.61))+0.369*exp(-0.12*(Y[2]+23.84)));
   dY[1] = (r_infinity-Y[1])/tau_r;
   i_CaL = Cm*g_CaL*(Y[2]-E_CaL)*Y[13]*Y[15]*Y[14];
   i_CaT = Cm*g_CaT*(Y[2]-E_CaT)*Y[16]*Y[17];
   i_f_Na = Cm*0.3833*g_if*(Y[2]-E_Na)*pow(Y[22], 2.0);
   i_f_K = Cm*0.6167*g_if*(Y[2]-E_K)*pow(Y[22], 2.0);
   i_f = i_f_Na+i_f_K;
   i_st = Cm*g_st*(Y[2]-E_st)*Y[23]*Y[24];
   i_Kr = Cm*g_Kr*(Y[2]-E_K)*(0.6*Y[18]+0.4*Y[19])*Y[20];
   i_Ks = Cm*g_Ks*(Y[2]-E_Ks)*pow(Y[21], 2.0);
   i_NaK = Cm*i_NaK_max/((1.0+pow(Km_Kp/Ko, 1.2))*(1.0+pow(Km_Nap/Nai, 1.3))*(1.0+exp(-(Y[2]-E_Na+120.0)/30.0)));
   k32 = exp(Qn*Y[2]/(2.0*RTOnF));
   di = 1.0+Y[11]/Kci*(1.0+exp(-Qci*Y[2]/RTOnF)+Nai/Kcni)+Nai/K1ni*(1.0+Nai/K2ni*(1.0+Nai/K3ni));
   k14 = Nai/K1ni*Nai/K2ni*(1.0+Nai/K3ni)*exp(Qn*Y[2]/(2.0*RTOnF))/di;
   k12 = Y[11]/Kci*exp(-Qci*Y[2]/RTOnF)/di;
   k41 = exp(-Qn*Y[2]/(2.0*RTOnF));
   x2 = k32*k43*(k14+k12)+k41*k12*(k34+k32);
   d_o = 1.0+Cao/Kco*(1.0+exp(Qco*Y[2]/RTOnF))+Nao/K1no*(1.0+Nao/K2no*(1.0+Nao/K3no));
   k21 = Cao/Kco*exp(Qco*Y[2]/RTOnF)/d_o;
   k23 = Nao/K1no*Nao/K2no*(1.0+Nao/K3no)*exp(-Qn*Y[2]/(2.0*RTOnF))/d_o;
   x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
   x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
   x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);
   i_NaCa = Cm*kNaCa*(x2*k21-x1*k12)/(x1+x2+x3+x4);
   i_b_Ca = Cm*g_b_Ca*(Y[2]-E_CaL);
   i_b_Na = Cm*g_b_Na*(Y[2]-E_Na);
   dY[2] = -(i_CaL+i_CaT+i_f+i_st+i_Kr+i_Ks+i_to+i_sus+i_NaK+i_NaCa+i_b_Ca+i_b_Na)/Cm;
   delta_fTC = kf_TC*Y[12]*(1.0-Y[6])-kb_TC*Y[6];
   dY[6] = delta_fTC;
   delta_fTMC = kf_TMC*Y[12]*(1.0-(Y[7]+Y[8]))-kb_TMC*Y[7];
   dY[7] = delta_fTMC;
   delta_fTMM = kf_TMM*Mgi*(1.0-(Y[7]+Y[8]))-kb_TMM*Y[8];
   dY[8] = delta_fTMM;
   delta_fCMi = kf_CM*Y[12]*(1.0-Y[3])-kb_CM*Y[3];
   dY[3] = delta_fCMi;
   delta_fCMs = kf_CM*Y[11]*(1.0-Y[4])-kb_CM*Y[4];
   dY[4] = delta_fCMs;
   delta_fCQ = kf_CQ*Y[9]*(1.0-Y[5])-kb_CQ*Y[5];
   dY[5] = delta_fCQ;
   j_Ca_dif = (Y[11]-Y[12])/tau_dif_Ca;
   j_up = P_up/(1.0+K_up/Y[12]);
   dY[12] = (j_Ca_dif*V_sub-j_up*V_nsr)/V_i-(CM_tot*delta_fCMi+TC_tot*delta_fTC+TMC_tot*delta_fTMC);
   j_SRCarel = ks*Y[26]*(Y[9]-Y[11]);
   dY[11] = j_SRCarel*V_jsr/V_sub-((i_CaL+i_CaT+i_b_Ca-2.0*i_NaCa)/(2.0*F*V_sub)+j_Ca_dif+CM_tot*delta_fCMs);
   j_tr = (Y[10]-Y[9])/tau_tr;
   dY[10] = j_up-j_tr*V_jsr/V_nsr;
   dY[9] = j_tr-(j_SRCarel+CQ_tot*delta_fCQ);
   dL_infinity = 1.0/(1.0+exp(-(Y[2]+13.5)/6.0));

   if (Y[2] == -35.0)
      adVm = -35.00001;
   else if (Y[2] == 0.0)
      adVm = 0.00001;
   else
      adVm = Y[2];

   alpha_dL = -0.02839*(adVm+35.0)/(exp(-(adVm+35.0)/2.5)-1.0)-0.0849*adVm/(exp(-adVm/4.8)-1.0);

   if (Y[2] == 5.0)
      bdVm = 5.00001;
   else
      bdVm = Y[2];

   beta_dL = 0.01143*(bdVm-5.0)/(exp((bdVm-5.0)/2.5)-1.0);
   tau_dL = 1.0/(alpha_dL+beta_dL);
   dY[13] = (dL_infinity-Y[13])/tau_dL;
   fCa_infinity = Km_fCa/(Km_fCa+Y[11]);
   tau_fCa = fCa_infinity/alpha_fCa;
   dY[14] = (fCa_infinity-Y[14])/tau_fCa;
   fL_infinity = 1.0/(1.0+exp((Y[2]+35.0)/7.3));
   tau_fL = 44.3+257.1*exp(-pow((Y[2]+32.5)/13.9, 2.0));
   dY[15] = (fL_infinity-Y[15])/tau_fL;
   dT_infinity = 1.0/(1.0+exp(-(Y[2]+26.3)/6.0));
   tau_dT = 1.0/(1.068*exp((Y[2]+26.3)/30.0)+1.068*exp(-(Y[2]+26.3)/30.0));
   dY[16] = (dT_infinity-Y[16])/tau_dT;
   fT_infinity = 1.0/(1.0+exp((Y[2]+61.7)/5.6));
   tau_fT = 1.0/(0.0153*exp(-(Y[2]+61.7)/83.3)+0.015*exp((Y[2]+61.7)/15.38));
   dY[17] = (fT_infinity-Y[17])/tau_fT;
   pa_infinity = 1.0/(1.0+exp(-(Y[2]+23.2)/10.6));
   tau_paS = 0.84655354/(0.0042*exp(Y[2]/17.0)+0.00015*exp(-Y[2]/21.6));
   tau_paF = 0.84655354/(0.0372*exp(Y[2]/15.9)+0.00096*exp(-Y[2]/22.5));
   dY[19] = (pa_infinity-Y[19])/tau_paS;
   dY[18] = (pa_infinity-Y[18])/tau_paF;
   pi_infinity = 1.0/(1.0+exp((Y[2]+28.6)/17.1));
   tau_pi = 1.0/(0.1*exp(-Y[2]/54.645)+0.656*exp(Y[2]/106.157));
   dY[20] = (pi_infinity-Y[20])/tau_pi;
   alpha_n = 0.014/(1.0+exp(-(Y[2]-40.0)/9.0));
   beta_n = 0.001*exp(-Y[2]/45.0);
   n_infinity = alpha_n/(alpha_n+beta_n);
   tau_n = 1.0/(alpha_n+beta_n);
   dY[21] = (n_infinity-Y[21])/tau_n;
   y_infinity = 1.0/(1.0+exp((Y[2]-VIf_half)/13.5));
   tau_y = 0.7166529/(exp(-(Y[2]+386.9)/45.302)+exp((Y[2]-73.08)/19.231));
   dY[22] = (y_infinity-Y[22])/tau_y;
   qa_infinity = 1.0/(1.0+exp(-(Y[2]+57.0)/5.0));
   alpha_qa = 1.0/(0.15*exp(-Y[2]/11.0)+0.2*exp(-Y[2]/700.0));
   beta_qa = 1.0/(16.0*exp(Y[2]/8.0)+15.0*exp(Y[2]/50.0));
   tau_qa = 1.0/(alpha_qa+beta_qa);
   dY[23] = (qa_infinity-Y[23])/tau_qa;
   alpha_qi = 1.0/(3100.0*exp(Y[2]/13.0)+700.0*exp(Y[2]/70.0));
   beta_qi = 1.0/(95.0*exp(-Y[2]/10.0)+50.0*exp(-Y[2]/700.0))+0.000229/(1.0+exp(-Y[2]/5.0));
   qi_infinity = alpha_qi/(alpha_qi+beta_qi);
   tau_qi = 6.65/(alpha_qi+beta_qi);
   dY[24] = (qi_infinity-Y[24])/tau_qi;
   kCaSR = MaxSR-(MaxSR-MinSR)/(1.0+pow(EC50_SR/Y[9], HSR));
   koSRCa = koCa/kCaSR;
   kiSRCa = kiCa*kCaSR;
   dY[27] = kim*Y[28]-kiSRCa*Y[11]*Y[27]-(koSRCa*pow(Y[11], 2.0)*Y[27]-kom*Y[26]);
   dY[26] = koSRCa*pow(Y[11], 2.0)*Y[27]-kom*Y[26]-(kiSRCa*Y[11]*Y[26]-kim*Y[25]);
   dY[25] = kiSRCa*Y[11]*Y[26]-kim*Y[25]-(kom*Y[25]-koSRCa*pow(Y[11], 2.0)*Y[28]);
   dY[28] = kom*Y[25]-koSRCa*pow(Y[11], 2.0)*Y[28]-(kim*Y[28]-kiSRCa*Y[11]*Y[27]);
}

//==============================================================================
// End of file
//==============================================================================
