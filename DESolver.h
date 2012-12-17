#ifndef DESOLVER_H
#define DESOLVER_H
#include <Cell.h>
#include <GnuplotPiper.h>
#include <stdio.h>
using namespace std;

/** @brief This is a general differential equations solver class
  */
class DESolver{
private:
    GnuplotPiper *Gpp;
    FILE *ofs;
public:
    DESolver(GnuplotPiper *GPP = NULL){
        Gpp = GPP;
        ofs = fopen("rst.m","w");
        fprintf(ofs,"clear all;\nx=[\n");
    }

    /** public constructor
        @param GPP - pointer to gnuplot piper, in vheart it is NULL by default
      */
    //DESolver(){
    //}

    /** integrate a cell over single "small" time step
        @param cell - pointer to the cell to be integrated
        @param dt - time step
        @param time - current time
      */
    virtual void IntegrateOverDt(Cell *cell, double dt, double time = 0, bool inSplitPhase = false) = 0;

    /** solves dynamical equations over given time period
        @param cell - pointer to the cell to be integrated
        @param MaxTime - time period in milliseconds
        @param dt - time step
      */
    void Solve(Cell *cell, double MaxTime = 1000, double dt = 0.005){
        unsigned long int MT = (unsigned long int)MaxTime/dt;
        for (unsigned long int i = 1; i <= MT; i++){
            double time = i*dt;
            unsigned long int step = MT/1000;
            this->IntegrateOverDt(cell,dt,time);
            if (Gpp) Gpp->pipe1D(i,cell->getV());

            if (i/step*step == i) {
               fprintf(ofs,"%g\t%g\t%g\n",time,cell->getV(),cell->getTension());
            }

        }
        fprintf(ofs,"];\nplot(x(:,1),x(:,2));\n");
        fprintf(ofs,"hold;\nplot(x(:,1),x(:,3),'r');\n");
        fprintf(stderr,"Solve finished\n");
    }
};

/** @brief class implementing forward Euler integration method
  */
class ForwardEulerSolver:public DESolver{

public:
    long count;
    //ForwardEulerSolver():DESolver(){}
    ForwardEulerSolver(GnuplotPiper *GPP = NULL):DESolver(GPP){count = 0;}
    /** integrate a cell over single "small" time step
        @param cell - pointer to the cell to be integrated
        @param dt - time step
        @param time - current time
      */
    void IntegrateOverDt(Cell *cell, double dt, double time = 0, bool inSplitPhase = false){
        cell->computeRightHandSide(time);

        double dV = cell->getDV();
        if (fabs(dV) > cell->splitVoltageDerivativeThreshold && !inSplitPhase){
            double split_dt = dt/cell->splitStepFactor;
            for (int i=0; i<int(cell->splitStepFactor); i++){
                this->IntegrateOverDt(cell,split_dt,time+split_dt*i,true);
            }
        }else{
            for (int i=0; i<cell->systemSize; i++){
                cell->Vars[i] += dt*cell->dVars[i];
                count++;
            }
        }

    }
};

#endif // DESOLVER_H
