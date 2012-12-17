
#ifndef __LUO_RUDY_I_MODEL_1991_HPP__
#define __LUO_RUDY_I_MODEL_1991_HPP__
#include <Cell.h>
//------------------------------------------------------------------------------

#define __LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_ 8+ENABLE_MECHANO_ELECTRIC_COUPLING

//------------------------------------------------------------------------------

#include <string>

//------------------------------------------------------------------------------

class luo_rudy_I_model_1991 : public Cell {
   public:
      //------------------------------------------------------------------------
      // State variables
      //------------------------------------------------------------------------

      double Y[__LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_];
      double dY[__LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_];
      // 0: h (dimensionless) (in fast_sodium_current_h_gate)
      // 1: j (dimensionless) (in fast_sodium_current_j_gate)
      // 2: m (dimensionless) (in fast_sodium_current_m_gate)
      // 3: Cai (millimolar) (in intracellular_calcium_concentration)
      // 4: V (millivolt) (in membrane)
      // 5: d (dimensionless) (in slow_inward_current_d_gate)
      // 6: f (dimensionless) (in slow_inward_current_f_gate)
      // 7: X (dimensionless) (in time_dependent_potassium_current_X_gate)

      std::string YNames[__LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_];
      std::string YUnits[__LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_];
      std::string YComponents[__LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_];

      //------------------------------------------------------------------------
      // Constants
      //------------------------------------------------------------------------

      double E_b;   // millivolt (in background_current)
      double g_b;   // milliS_per_cm2 (in background_current)
      double g_Na;   // milliS_per_cm2 (in fast_sodium_current)
      double Ki;   // millimolar (in ionic_concentrations)
      double Ko;   // millimolar (in ionic_concentrations)
      double Nai;   // millimolar (in ionic_concentrations)
      double Nao;   // millimolar (in ionic_concentrations)
      double C;   // microF_per_cm2 (in membrane)
      double F;   // coulomb_per_mole (in membrane)
      double R;   // joule_per_kilomole_kelvin (in membrane)
      double T;   // kelvin (in membrane)
      double stim_amplitude;   // microA_per_cm2 (in membrane)
      double stim_duration;   // millisecond (in membrane)
      double stim_end;   // millisecond (in membrane)
      double stim_period;   // millisecond (in membrane)
      double stim_start;   // millisecond (in membrane)
      double g_Kp;   // milliS_per_cm2 (in plateau_potassium_current)
      double PR_NaK;   // dimensionless (in time_dependent_potassium_current)

      //------------------------------------------------------------------------
      // Computed variables
      //------------------------------------------------------------------------

      double i_b;   // microA_per_cm2 (in background_current)
      double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
      double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
      double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
      double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
      double alpha_m;   // per_millisecond (in fast_sodium_current_m_gate)
      double beta_m;   // per_millisecond (in fast_sodium_current_m_gate)
      double E_Na;   // millivolt (in fast_sodium_current)
      double i_Na;   // microA_per_cm2 (in fast_sodium_current)
      double I_stim;   // microA_per_cm2 (in membrane)
      double E_Kp;   // millivolt (in plateau_potassium_current)
      double Kp;   // dimensionless (in plateau_potassium_current)
      double i_Kp;   // microA_per_cm2 (in plateau_potassium_current)
      double alpha_d;   // per_millisecond (in slow_inward_current_d_gate)
      double beta_d;   // per_millisecond (in slow_inward_current_d_gate)
      double alpha_f;   // per_millisecond (in slow_inward_current_f_gate)
      double beta_f;   // per_millisecond (in slow_inward_current_f_gate)
      double E_si;   // millivolt (in slow_inward_current)
      double i_si;   // microA_per_cm2 (in slow_inward_current)
      double alpha_X;   // per_millisecond (in time_dependent_potassium_current_X_gate)
      double beta_X;   // per_millisecond (in time_dependent_potassium_current_X_gate)
      double Xi;   // dimensionless (in time_dependent_potassium_current_Xi_gate)
      double E_K;   // millivolt (in time_dependent_potassium_current)
      double g_K;   // milliS_per_cm2 (in time_dependent_potassium_current)
      double i_K;   // microA_per_cm2 (in time_dependent_potassium_current)
      double K1_infinity;   // dimensionless (in time_independent_potassium_current_K1_gate)
      double alpha_K1;   // per_millisecond (in time_independent_potassium_current_K1_gate)
      double beta_K1;   // per_millisecond (in time_independent_potassium_current_K1_gate)
      double E_K1;   // millivolt (in time_independent_potassium_current)
      double g_K1;   // milliS_per_cm2 (in time_independent_potassium_current)
      double i_K1;   // microA_per_cm2 (in time_independent_potassium_current)
      double g_Si;
      //------------------------------------------------------------------------
      // Procedures
      //------------------------------------------------------------------------

      void init();
      void compute(double time);

      double getV(){
          return this->Y[4];
      }

      double getDV(){
          return this->dY[4];
      }

      double getStim(double time){
          if ((time-floor(time/stim_period)*stim_period >= stim_start) &&
                  (time-floor(time/stim_period)*stim_period <= stim_start+stim_duration))
              return -stim_amplitude;
          else
              return 0.0;
      }

      void incrementV(double value){
          this->Y[4] += value;
      }

      void setGK1(double value){
          g_K1 = value*sqrt(Ko/5.4);
      }

      void setGsi(double value){
          g_Si = value;
      }

      void setGkp(double value){
          g_Kp = value;
      }
      void setGk(double value){
        g_K = value*sqrt(Ko/5.4);
      }
      //------------------------------------------------------------------
      // Public Constructor
      //------------------------------------------------------------------
      luo_rudy_I_model_1991(bool enableMech = false) : Cell(this->Y, this->dY, __LUO_RUDY_I_MODEL_1991_NB_OF_STATE_VARIABLES_,
                                     &this->stim_amplitude, &this->stim_duration, &this->stim_end,
                                     &this->stim_period, &this->stim_start, enableMech){}

};

//------------------------------------------------------------------------------

#endif

//==============================================================================
// End of file
//==============================================================================
