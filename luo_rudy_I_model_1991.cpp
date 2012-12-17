//==============================================================================
// CellML file:   C:\Program Files\COR\Models\luo_rudy_I_model_1991.cellml
// CellML model:  luo_rudy_I_model_1991
// Date and time: 23.03.2012 at 16:12:20
//------------------------------------------------------------------------------
// Conversion from CellML 1.0 to C++ was done using COR (0.9.31.1409)
//    Copyright 2002-2012 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================

#include "luo_rudy_I_model_1991.h"

//------------------------------------------------------------------------------

#include <math.h>

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void luo_rudy_I_model_1991::init()
{
   //---------------------------------------------------------------------------
   // State variables
   //---------------------------------------------------------------------------

   Y[0] = 0.982660523699656;   // h (dimensionless) (in fast_sodium_current_h_gate)
   Y[1] = 0.989108212766685;   // j (dimensionless) (in fast_sodium_current_j_gate)
   Y[2] = 0.00171338077730188;   // m (dimensionless) (in fast_sodium_current_m_gate)
   Y[3] = 0.00017948816388306;   // Cai (millimolar) (in intracellular_calcium_concentration)
   Y[4] = -84.3801107371;   // V (millivolt) (in membrane)
   Y[5] = 0.00302126301779861;   // d (dimensionless) (in slow_inward_current_d_gate)
   Y[6] = 0.999967936476325;   // f (dimensionless) (in slow_inward_current_f_gate)
   Y[7] = 0.0417603108167287;   // X (dimensionless) (in time_dependent_potassium_current_X_gate)

   YNames[0].assign("h");
   YNames[1].assign("j");
   YNames[2].assign("m");
   YNames[3].assign("Cai");
   YNames[4].assign("V");
   YNames[5].assign("d");
   YNames[6].assign("f");
   YNames[7].assign("X");

   YUnits[0].assign("dimensionless");
   YUnits[1].assign("dimensionless");
   YUnits[2].assign("dimensionless");
   YUnits[3].assign("millimolar");
   YUnits[4].assign("millivolt");
   YUnits[5].assign("dimensionless");
   YUnits[6].assign("dimensionless");
   YUnits[7].assign("dimensionless");

   YComponents[0].assign("fast_sodium_current_h_gate");
   YComponents[1].assign("fast_sodium_current_j_gate");
   YComponents[2].assign("fast_sodium_current_m_gate");
   YComponents[3].assign("intracellular_calcium_concentration");
   YComponents[4].assign("membrane");
   YComponents[5].assign("slow_inward_current_d_gate");
   YComponents[6].assign("slow_inward_current_f_gate");
   YComponents[7].assign("time_dependent_potassium_current_X_gate");

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   E_b = -59.87;   // millivolt (in background_current)
   g_b = 0.03921;   // milliS_per_cm2 (in background_current)
   g_Na = 23.0;   // milliS_per_cm2 (in fast_sodium_current)
   Ki = 145.0;   // millimolar (in ionic_concentrations)
   Ko = 5.4;   // millimolar (in ionic_concentrations)
   Nai = 18.0;   // millimolar (in ionic_concentrations)
   Nao = 140.0;   // millimolar (in ionic_concentrations)
   C = 1.0;   // microF_per_cm2 (in membrane)
   F = 96484.6;   // coulomb_per_mole (in membrane)
   R = 8314.0;   // joule_per_kilomole_kelvin (in membrane)
   T = 310.0;   // kelvin (in membrane)
   stim_amplitude = -25.5;   // microA_per_cm2 (in membrane)
   stim_duration = 2.0;   // millisecond (in membrane)
   stim_end = 9000.0;   // millisecond (in membrane)
   stim_period = 1000.0;   // millisecond (in membrane)
   stim_start = 100.0;   // millisecond (in membrane)
   g_Kp = 0.0183;   // milliS_per_cm2 (in plateau_potassium_current)
   PR_NaK = 0.01833;   // dimensionless (in time_dependent_potassium_current)

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   E_Na = R*T/F*log(Nao/Nai);
   g_K = 0.282*sqrt(Ko/5.4);
   E_K = R*T/F*log((Ko+PR_NaK*Nao)/(Ki+PR_NaK*Nai));
   g_K1 = 0.6047*sqrt(Ko/5.4);
   E_K1 = R*T/F*log(Ko/Ki);
   E_Kp = E_K1;
   g_Si = 0.09;
   //--------------------------------------------------------------------------
   // Integration params
   //--------------------------------------------------------------------------
   this->dt = 0.03;
   this->splitStepFactor = 10; //this implies minimal step of 0.001
   this->splitVoltageDerivativeThreshold = 2;
}

//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

void luo_rudy_I_model_1991::compute(double time)
{
   // time: time (millisecond)

   i_b = g_b*(Y[4]-E_b);
   i_Na = g_Na*pow(Y[2], 3.0)*Y[0]*Y[1]*(Y[4]-E_Na);

   if (Y[4] < -40.0)
      alpha_h = 0.135*exp((80.0+Y[4])/-6.8);
   else
      alpha_h = 0.0;

   if (Y[4] < -40.0)
      beta_h = 3.56*exp(0.079*Y[4])+310000.0*exp(0.35*Y[4]);
   else
      beta_h = 1.0/(0.13*(1.0+exp((Y[4]+10.66)/-11.1)));

   dY[0] = alpha_h*(1.0-Y[0])-beta_h*Y[0];

   if (Y[4] < -40.0)
      alpha_j = (-127140.0*exp(0.2444*Y[4])-0.00003474*exp(-0.04391*Y[4]))*(Y[4]+37.78)/(1.0+exp(0.311*(Y[4]+79.23)));
   else
      alpha_j = 0.0;

   if (Y[4] < -40.0)
      beta_j = 0.1212*exp(-0.01052*Y[4])/(1.0+exp(-0.1378*(Y[4]+40.14)));
   else
      beta_j = 0.3*exp(-0.0000002535*Y[4])/(1.0+exp(-0.1*(Y[4]+32.0)));

   dY[1] = alpha_j*(1.0-Y[1])-beta_j*Y[1];
   alpha_m = 0.32*(Y[4]+47.13)/(1.0-exp(-0.1*(Y[4]+47.13)));
   beta_m = 0.08*exp(-Y[4]/11.0);
   dY[2] = alpha_m*(1.0-Y[2])-beta_m*Y[2];
   E_si = 7.7-13.0287*log(Y[3]/1.0);
   i_si = g_Si*Y[5]*Y[6]*(Y[4]-E_si);
   dY[3] = -0.0001/1.0*i_si+0.07*(0.0001-Y[3]);

   if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
      I_stim = stim_amplitude;
   else
      I_stim = 0.0;

   if (Y[4] > -100.0)
      Xi = 2.837*(exp(0.04*(Y[4]+77.0))-1.0)/((Y[4]+77.0)*exp(0.04*(Y[4]+35.0)));
   else
      Xi = 1.0;

   i_K = g_K*Y[7]*Xi*(Y[4]-E_K);
   alpha_K1 = 1.02/(1.0+exp(0.2385*(Y[4]-E_K1-59.215)));
   beta_K1 = (0.49124*exp(0.08032*(Y[4]+5.476-E_K1))+1.0*exp(0.06175*(Y[4]-(E_K1+594.31))))/(1.0+exp(-0.5143*(Y[4]-E_K1+4.753)));
   K1_infinity = alpha_K1/(alpha_K1+beta_K1);
   i_K1 = g_K1*K1_infinity*(Y[4]-E_K1);
   Kp = 1.0/(1.0+exp((7.488-Y[4])/5.98));
   i_Kp = g_Kp*Kp*(Y[4]-E_Kp);
   dY[4] = -1.0/C*(I_stim+i_Na+i_si+i_K+i_K1+i_Kp+i_b+calcISAC());
   alpha_d = 0.095*exp(-0.01*(Y[4]-5.0))/(1.0+exp(-0.072*(Y[4]-5.0)));
   beta_d = 0.07*exp(-0.017*(Y[4]+44.0))/(1.0+exp(0.05*(Y[4]+44.0)));
   dY[5] = alpha_d*(1.0-Y[5])-beta_d*Y[5];
   alpha_f = 0.012*exp(-0.008*(Y[4]+28.0))/(1.0+exp(0.15*(Y[4]+28.0)));
   beta_f = 0.0065*exp(-0.02*(Y[4]+30.0))/(1.0+exp(-0.2*(Y[4]+30.0)));
   dY[6] = alpha_f*(1.0-Y[6])-beta_f*Y[6];
   alpha_X = 0.0005*exp(0.083*(Y[4]+50.0))/(1.0+exp(0.057*(Y[4]+50.0)));
   beta_X = 0.0013*exp(-0.06*(Y[4]+20.0))/(1.0+exp(-0.04*(Y[4]+20.0)));
   dY[7] = alpha_X*(1.0-Y[7])-beta_X*Y[7];
}

//==============================================================================
// End of file
//==============================================================================
