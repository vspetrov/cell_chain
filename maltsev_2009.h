
#ifndef __MALTSEV_2009_HPP__
#define __MALTSEV_2009_HPP__
#include <Cell.h>
#include <stdio.h>
#include <stdlib.h>
//------------------------------------------------------------------------------

#define _MALTSEV_2009_NB_OF_STATE_VARIABLES_ 29

//------------------------------------------------------------------------------

#include <string>

//------------------------------------------------------------------------------

class maltsev_2009 : public Cell
{
   public:
      //------------------------------------------------------------------------
      // State variables
      //------------------------------------------------------------------------

      double Y[_MALTSEV_2009_NB_OF_STATE_VARIABLES_];
      double dY[_MALTSEV_2009_NB_OF_STATE_VARIABLES_];
      // 0: q (dimensionless) (in AP_sensitive_currents_q_gate)
      // 1: r (dimensionless) (in AP_sensitive_currents_r_gate)
      // 2: Vm (millivolt) (in Vm)
      // 3: fCMi (dimensionless) (in calcium_buffering)
      // 4: fCMs (dimensionless) (in calcium_buffering)
      // 5: fCQ (dimensionless) (in calcium_buffering)
      // 6: fTC (dimensionless) (in calcium_buffering)
      // 7: fTMC (dimensionless) (in calcium_buffering)
      // 8: fTMM (dimensionless) (in calcium_buffering)
      // 9: Ca_jsr (millimolar) (in calcium_dynamics)
      // 10: Ca_nsr (millimolar) (in calcium_dynamics)
      // 11: Ca_sub (millimolar) (in calcium_dynamics)
      // 12: Cai (millimolar) (in calcium_dynamics)
      // 13: dL (dimensionless) (in i_CaL_dL_gate)
      // 14: fCa (dimensionless) (in i_CaL_fCa_gate)
      // 15: fL (dimensionless) (in i_CaL_fL_gate)
      // 16: dT (dimensionless) (in i_CaT_dT_gate)
      // 17: fT (dimensionless) (in i_CaT_fT_gate)
      // 18: paF (dimensionless) (in i_Kr_pa_gate)
      // 19: paS (dimensionless) (in i_Kr_pa_gate)
      // 20: pi_ (dimensionless) (in i_Kr_pi_gate)
      // 21: n (dimensionless) (in i_Ks_n_gate)
      // 22: y (dimensionless) (in i_f_y_gate)
      // 23: qa (dimensionless) (in i_st_qa_gate)
      // 24: qi (dimensionless) (in i_st_qi_gate)
      // 25: I (dimensionless) (in j_SRCarel)
      // 26: O (dimensionless) (in j_SRCarel)
      // 27: R1 (dimensionless) (R in j_SRCarel)
      // 28: RI (dimensionless) (in j_SRCarel)

      std::string YNames[_MALTSEV_2009_NB_OF_STATE_VARIABLES_];
      std::string YUnits[_MALTSEV_2009_NB_OF_STATE_VARIABLES_];
      std::string YComponents[_MALTSEV_2009_NB_OF_STATE_VARIABLES_];

      //------------------------------------------------------------------------
      // Constants
      //------------------------------------------------------------------------

      double g_sus;   // nanoS_per_picoF (in AP_sensitive_currents)
      double g_to;   // nanoS_per_picoF (in AP_sensitive_currents)
      double CM_tot;   // millimolar (in calcium_buffering)
      double CQ_tot;   // millimolar (in calcium_buffering)
      double TC_tot;   // millimolar (in calcium_buffering)
      double TMC_tot;   // millimolar (in calcium_buffering)
      double kb_CM;   // per_millisecond (in calcium_buffering)
      double kb_CQ;   // per_millisecond (in calcium_buffering)
      double kb_TC;   // per_millisecond (in calcium_buffering)
      double kb_TMC;   // per_millisecond (in calcium_buffering)
      double kb_TMM;   // per_millisecond (in calcium_buffering)
      double kf_CM;   // per_millimolar_millisecond (in calcium_buffering)
      double kf_CQ;   // per_millimolar_millisecond (in calcium_buffering)
      double kf_TC;   // per_millimolar_millisecond (in calcium_buffering)
      double kf_TMC;   // per_millimolar_millisecond (in calcium_buffering)
      double kf_TMM;   // per_millimolar_millisecond (in calcium_buffering)
      double Km_fCa;   // millimolar (in i_CaL_fCa_gate)
      double alpha_fCa;   // per_millisecond (in i_CaL_fCa_gate)
      double E_CaL;   // millivolt (in i_CaL)
      double g_CaL;   // nanoS_per_picoF (in i_CaL)
      double E_CaT;   // millivolt (in i_CaT)
      double g_CaT;   // nanoS_per_picoF (in i_CaT)
      double g_Kr;   // nanoS_per_picoF (in i_Kr)
      double g_Ks;   // nanoS_per_picoF (in i_Ks)
      double K1ni;   // millimolar (in i_NaCa)
      double K1no;   // millimolar (in i_NaCa)
      double K2ni;   // millimolar (in i_NaCa)
      double K2no;   // millimolar (in i_NaCa)
      double K3ni;   // millimolar (in i_NaCa)
      double K3no;   // millimolar (in i_NaCa)
      double Kci;   // millimolar (in i_NaCa)
      double Kcni;   // millimolar (in i_NaCa)
      double Kco;   // millimolar (in i_NaCa)
      double Qci;   // dimensionless (in i_NaCa)
      double Qco;   // dimensionless (in i_NaCa)
      double Qn;   // dimensionless (in i_NaCa)
      double kNaCa;   // picoA_per_picoF (in i_NaCa)
      double Km_Kp;   // millimolar (in i_NaK)
      double Km_Nap;   // millimolar (in i_NaK)
      double i_NaK_max;   // picoA_per_picoF (in i_NaK)
      double g_b_Ca;   // nanoS_per_picoF (in i_b_Ca)
      double g_b_Na;   // nanoS_per_picoF (in i_b_Na)
      double VIf_half;   // millivolt (in i_f_y_gate)
      double g_if;   // nanoS_per_picoF (in i_f)
      double E_st;   // millivolt (in i_st)
      double g_st;   // nanoS_per_picoF (in i_st)
      double K_up;   // millimolar (in intracellular_calcium_fluxes)
      double P_up;   // millimolar_per_millisecond (in intracellular_calcium_fluxes)
      double tau_dif_Ca;   // millisecond (in intracellular_calcium_fluxes)
      double tau_tr;   // millisecond (in intracellular_calcium_fluxes)
      double EC50_SR;   // millimolar (in j_SRCarel)
      double HSR;   // dimensionless (in j_SRCarel)
      double MaxSR;   // dimensionless (in j_SRCarel)
      double MinSR;   // dimensionless (in j_SRCarel)
      double kiCa;   // per_millimolar_millisecond (in j_SRCarel)
      double kim;   // per_millisecond (in j_SRCarel)
      double koCa;   // per_millimolar2_millisecond (in j_SRCarel)
      double kom;   // per_millisecond (in j_SRCarel)
      double ks;   // per_millisecond (in j_SRCarel)
      double Cao;   // millimolar (in model_parameters)
      double Cm;   // picoF (in model_parameters)
      double F;   // coulomb_per_mole (in model_parameters)
      double Ki;   // millimolar (in model_parameters)
      double Ko;   // millimolar (in model_parameters)
      double L_cell;   // micrometre (in model_parameters)
      double L_sub;   // micrometre (in model_parameters)
      double Mgi;   // millimolar (in model_parameters)
      double Nai;   // millimolar (in model_parameters)
      double Nao;   // millimolar (in model_parameters)
      double R2;   // joule_per_kilomole_kelvin (R in model_parameters)
      double R_cell;   // micrometre (in model_parameters)
      double T;   // kelvin (in model_parameters)
      double V_i_part;   // dimensionless (in model_parameters)
      double V_jsr_part;   // dimensionless (in model_parameters)
      double V_nsr_part;   // dimensionless (in model_parameters)
      double stim_amplitude;   // microA_per_cm2 (in membrane)
      double stim_duration;   // millisecond (in membrane)
      double stim_end;   // millisecond (in membrane)
      double stim_period;   // millisecond (in membrane)
      double stim_start;   // millisecond (in membrane)
      //------------------------------------------------------------------------
      // Computed variables
      //------------------------------------------------------------------------

      double q_infinity;   // dimensionless (in AP_sensitive_currents_q_gate)
      double tau_q;   // millisecond (in AP_sensitive_currents_q_gate)
      double r_infinity;   // dimensionless (in AP_sensitive_currents_r_gate)
      double tau_r;   // millisecond (in AP_sensitive_currents_r_gate)
      double i_sus;   // picoA (in AP_sensitive_currents)
      double i_to;   // picoA (in AP_sensitive_currents)
      double delta_fCMi;   // per_millisecond (in calcium_buffering)
      double delta_fCMs;   // per_millisecond (in calcium_buffering)
      double delta_fCQ;   // per_millisecond (in calcium_buffering)
      double delta_fTC;   // per_millisecond (in calcium_buffering)
      double delta_fTMC;   // per_millisecond (in calcium_buffering)
      double delta_fTMM;   // per_millisecond (in calcium_buffering)
      double E_K;   // millivolt (in electric_potentials)
      double E_Ks;   // millivolt (in electric_potentials)
      double E_Na;   // millivolt (in electric_potentials)
      double adVm;   // millivolt (in i_CaL_dL_gate)
      double alpha_dL;   // per_millisecond (in i_CaL_dL_gate)
      double bdVm;   // millivolt (in i_CaL_dL_gate)
      double beta_dL;   // per_millisecond (in i_CaL_dL_gate)
      double dL_infinity;   // dimensionless (in i_CaL_dL_gate)
      double tau_dL;   // millisecond (in i_CaL_dL_gate)
      double fCa_infinity;   // dimensionless (in i_CaL_fCa_gate)
      double tau_fCa;   // millisecond (in i_CaL_fCa_gate)
      double fL_infinity;   // dimensionless (in i_CaL_fL_gate)
      double tau_fL;   // millisecond (in i_CaL_fL_gate)
      double i_CaL;   // picoA (in i_CaL)
      double dT_infinity;   // dimensionless (in i_CaT_dT_gate)
      double tau_dT;   // millisecond (in i_CaT_dT_gate)
      double fT_infinity;   // dimensionless (in i_CaT_fT_gate)
      double tau_fT;   // millisecond (in i_CaT_fT_gate)
      double i_CaT;   // picoA (in i_CaT)
      double pa_infinity;   // dimensionless (in i_Kr_pa_gate)
      double tau_paF;   // millisecond (in i_Kr_pa_gate)
      double tau_paS;   // millisecond (in i_Kr_pa_gate)
      double pi_infinity;   // dimensionless (in i_Kr_pi_gate)
      double tau_pi;   // millisecond (in i_Kr_pi_gate)
      double i_Kr;   // picoA (in i_Kr)
      double alpha_n;   // per_millisecond (in i_Ks_n_gate)
      double beta_n;   // per_millisecond (in i_Ks_n_gate)
      double n_infinity;   // dimensionless (in i_Ks_n_gate)
      double tau_n;   // millisecond (in i_Ks_n_gate)
      double i_Ks;   // picoA (in i_Ks)
      double RTOnF;   // millivolt (in i_NaCa)
      double di;   // dimensionless (in i_NaCa)
      double d_o;   // dimensionless (in i_NaCa)
      double i_NaCa;   // picoA (in i_NaCa)
      double k12;   // dimensionless (in i_NaCa)
      double k14;   // dimensionless (in i_NaCa)
      double k21;   // dimensionless (in i_NaCa)
      double k23;   // dimensionless (in i_NaCa)
      double k32;   // dimensionless (in i_NaCa)
      double k34;   // dimensionless (in i_NaCa)
      double k41;   // dimensionless (in i_NaCa)
      double k43;   // dimensionless (in i_NaCa)
      double x1;   // dimensionless (in i_NaCa)
      double x2;   // dimensionless (in i_NaCa)
      double x3;   // dimensionless (in i_NaCa)
      double x4;   // dimensionless (in i_NaCa)
      double i_NaK;   // picoA (in i_NaK)
      double i_b_Ca;   // picoA (in i_b_Ca)
      double i_b_Na;   // picoA (in i_b_Na)
      double tau_y;   // millisecond (in i_f_y_gate)
      double y_infinity;   // dimensionless (in i_f_y_gate)
      double i_f;   // picoA (in i_f)
      double i_f_K;   // picoA (in i_f)
      double i_f_Na;   // picoA (in i_f)
      double alpha_qa;   // per_millisecond (in i_st_qa_gate)
      double beta_qa;   // per_millisecond (in i_st_qa_gate)
      double qa_infinity;   // dimensionless (in i_st_qa_gate)
      double tau_qa;   // millisecond (in i_st_qa_gate)
      double alpha_qi;   // per_millisecond (in i_st_qi_gate)
      double beta_qi;   // per_millisecond (in i_st_qi_gate)
      double qi_infinity;   // dimensionless (in i_st_qi_gate)
      double tau_qi;   // millisecond (in i_st_qi_gate)
      double i_st;   // picoA (in i_st)
      double j_Ca_dif;   // millimolar_per_millisecond (in intracellular_calcium_fluxes)
      double j_tr;   // millimolar_per_millisecond (in intracellular_calcium_fluxes)
      double j_up;   // millimolar_per_millisecond (in intracellular_calcium_fluxes)
      double j_SRCarel;   // millimolar_per_millisecond (in j_SRCarel)
      double kCaSR;   // dimensionless (in j_SRCarel)
      double kiSRCa;   // per_millimolar_millisecond (in j_SRCarel)
      double koSRCa;   // per_millimolar2_millisecond (in j_SRCarel)
      double V_cell;   // picolitre (in model_parameters)
      double V_i;   // picolitre (in model_parameters)
      double V_jsr;   // picolitre (in model_parameters)
      double V_nsr;   // picolitre (in model_parameters)
      double V_sub;   // picolitre (in model_parameters)

      //------------------------------------------------------------------------
      // Procedures
      //------------------------------------------------------------------------

      void init();
      void compute(double time);

      double getV(){
          return this->Y[2];
      }

      double getDV(){
          return this->dY[2];
      }

      double getStim(double time){
          if ((time-floor(time/stim_period)*stim_period >= stim_start) &&
              (time-floor(time/stim_period)*stim_period <= stim_start+stim_duration))
              return -stim_amplitude;
          else
              return 0.0;
      }

      void incrementV(double value){
          this->Y[2] += value;
      }
      //------------------------------------------------------------------
      // Public Constructor
      //------------------------------------------------------------------
      maltsev_2009() : Cell(this->Y, this->dY, _MALTSEV_2009_NB_OF_STATE_VARIABLES_,
                            &this->stim_amplitude, &this->stim_duration, &this->stim_end,
                            &this->stim_period, &this->stim_start){
          //--------------------------------------------------------------------------
          // Integration params
          //--------------------------------------------------------------------------
          this->dt = 0.01;
          this->splitStepFactor = 10; //this implies minimal step of 0.001
          this->splitVoltageDerivativeThreshold = 2;
      }

};

//------------------------------------------------------------------------------

#endif

//==============================================================================
// End of file
//==============================================================================
