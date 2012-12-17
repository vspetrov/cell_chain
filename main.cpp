#include "DESolver.h"
#include "luo_rudy_I_model_1991.h"
#include "maltsev_2009.h"
#include <cellchain.h>
#include <stdio.h>

static void integrateOneCell(double MaxTime){
    Cell *cell;
    cell = new luo_rudy_I_model_1991(true);
    cell->init();
    //cell->incrementV(50);
    *cell->stimulus.stim_amplitude = -20;

    GnuplotPiper Gpp(10,MaxTime,cell->dt,100);
    ForwardEulerSolver Solver(&Gpp);
    Solver.Solve(cell,MaxTime,cell->dt);
    printf("Number of steps: %li\n",Solver.count);
    delete cell;
}

int main(int argc, char *argv[])
{

    integrateOneCell(5000);
//    vector<int> toPace;
//    toPace.push_back(1);
//    toPace.push_back(2);
//    toPace.push_back(3);
//    CellChain chain(205,15,30,20,110,30);
//    chain.setStimStartColl(10);
//    chain.setStimAmplitudeColl(0);
//   // chain.setStimAmplitudeColl(-25,toPace);
//    chain.Solve(40000, 30000);

    fprintf(stderr,"Solve finished\n");
    return 0;
}
