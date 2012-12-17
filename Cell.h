#ifndef CELL_H
#define CELL_H
#include <math.h>
#include <string>
#include <iostream>
/** @brief This structure describes the external stimulus to the given Cell
 */
typedef struct stimulus{
    double *stim_amplitude;   /**< microA_per_cm2 (in membrane) */
    double *stim_duration;   /**< millisecond (in membrane) */
    double *stim_end;   /**< millisecond (in membrane) */
    double *stim_period;   /**< millisecond (in membrane) */
    double *stim_start;   /**< millisecond (in membrane) */
}stimulus_t;


#define ENABLE_MECHANO_ELECTRIC_COUPLING 1

/** @brief This is the main class for the Cell object. Objects of this type are used in virtualheart class.
    Each Electrical Model class must inherit from this Class.
 */
class Cell{
public:
    int systemSize; /**< the dimension of the dynamical system */
    double *Vars; /**< pointer to the array of dynamical variables */
    double *dVars; /**< pointer to the array of derivatives of the dynamical variables */
    stimulus_t stimulus; /**< stimulus to the Cell object */
    double tension;
    double *s_tension;
    double dt; /**integration step for this model*/
    double splitStepFactor; /**factor for integration algorithm with bariable time step*/
    double splitVoltageDerivativeThreshold; /**threshold when time step split is enabled*/
    /**  public constructor
         @param variables - pointer to the array of dynamical variables
         @param derivatives - pointer to the array of derivatives of the dynamical variables
         @param systemSize - the dimension of the dynamical system
         @param s_a - stimulus amplitude
         @param s_d - stimulus duration
         @param s_e - stimulus end time
         @param s_p - stimulus period
         @param s_s - stimulus start time
     */
    bool withMechanics;
    Cell(double *variables, double *derivatives, int systemSize,
         double *s_a, double *s_d, double *s_e, double *s_p, double *s_s, bool enableMech = false){
        this->Vars = variables;
        this->dVars = derivatives;
        this->systemSize = systemSize;
        std::cout << "System Size = " << systemSize << std::endl;
        stimulus.stim_amplitude = s_a;
        stimulus.stim_duration = s_d;
        stimulus.stim_end = s_e;
        stimulus.stim_period = s_p;
        stimulus.stim_start = s_s;
        withMechanics = enableMech;
#if ENABLE_MECHANO_ELECTRIC_COUPLING
        s_tension = &Vars[systemSize-1];
        *s_tension = 0; // Initialize tension variable to zero
        tension = 0;
#endif
    }

    /** initialize the cell: allocate memory, set up constants and initial conditions
     */
    virtual void init() = 0;

    /** calculate derivatives of the dynamical variables at the time moment t
        @param t - time moment
     */
    virtual void compute(double t) = 0;

    /** get the current voltage value of the cell
        @return  the voltage
     */
    virtual double getV() = 0;

    /** get the current voltage derivative value of the cell
        @return  the voltage
     */
    virtual double getDV() = 0;

    /** get the current voltage derivative value of the cell
        @return  the voltage
     */


    /** get stimulus amplitude value (depricated)
        @param time - time moment
        @return stimulus amplitude
      */
    virtual double getStim(double time) = 0;

    /** increment the value of the cells voltage by given value
        @param value - how much to increment
      */
    virtual void incrementV(double value) = 0;

    /** get stimulus info
        @return pointer to the stimulus structure that holds all stim info
      */
    virtual stimulus_t* stimAmplitude()
    {
        return &stimulus;
    }

    /** Support for Mechanics
      */
    static const double G_s_mech = 0.005;
    static const double E_mech_threshold = 20;
    static const double V_t = 10;
    static const double varepsilon_1 = 0.1;
    static const double varepsilon_2 = 0.01;
    static const double V_o = -80;
    static const double S_o = 0;

    virtual double calcISAC(){
        if (!withMechanics)
            return 0;

        tension = S_o - *s_tension;
        return G_s_mech * sigmoid(tension)*tension*(getV() - E_mech_threshold);
    }

    inline double sigmoid(double x){
        return 1./(1+exp(-x*10));
    }

    void computeMechanics(double V){
        double epsilon = varepsilon_1 + (varepsilon_2 - varepsilon_1)*sigmoid(getV()-V_o);
        dVars[systemSize-1] = epsilon*(getV()/V_t-*s_tension);
    }

    void computeRightHandSide(double time){
        compute(time);
        if (withMechanics){
            computeMechanics(this->getV());
        }
    }
    double getTension(){
        return tension;
    }
};

#endif // CELL_H
