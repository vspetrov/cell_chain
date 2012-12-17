#ifndef CELLCHAIN_H
#define CELLCHAIN_H

#include <Cell.h>
#include <DESolver.h>
#include <luo_rudy_I_model_1991.h>
#include <vector>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

class CellChain
{
private:
    int size;
    int sanLen;
    int atriaLen;
    int avLen;
    int purkLen;
    int ventLen;

    double *D;

    tbb::task_scheduler_init *tbb_init;
    void SolveDt(double dt, double time);
    void initD();

public:
    Cell **cells;
    DESolver *Solver;
    void Coupling(int i);
    CellChain(int size, int sanLen, int atriaLen, int avLen, int purkLen, int ventLen);
    Cell *getCell(int i);
    void Solve(double MaxTime, double skipTime);
    void setStimAmplitudeColl(double A, std::vector<int> ids = std::vector<int>());
    void setStimStartColl(double time);

};

class ParallelSolver{
    CellChain *chain;
    double dt;
    double time;
    int flag;//0 - IntegrateDT; 1  -Coupling
public:
    void operator()( const tbb::blocked_range<size_t>& r ) const {
        if (!flag)
            for( size_t i=r.begin(); i!=r.end(); ++i ){
                chain->Solver->IntegrateOverDt(chain->cells[i],chain->cells[i]->dt,time);
            }
        else
            for( size_t i=r.begin(); i!=r.end(); ++i ){
                chain->Coupling(i);
            }
    }
    ParallelSolver(CellChain *chain, double dt, double time, int flag):
        chain(chain),dt(dt),time(time),flag(flag){}
};

#endif // CELLCHAIN_H
