#include "cellchain.h"

CellChain::CellChain(int size, int sanLen, int atriaLen, int avLen, int purkLen, int ventLen)
{
    this->size = size;
    this->sanLen = sanLen;
    this->atriaLen = atriaLen;
    this->avLen = avLen;
    this->purkLen = purkLen;
    this->ventLen = ventLen;
    this->cells = new Cell*[size];
    this->D = new double[size];
    this->Solver = new ForwardEulerSolver(NULL);
    //init
    for (int i=0; i<size; i++){
        cells[i] = (Cell *)new luo_rudy_I_model_1991();
        cells[i]->init();
    }
    initD();

    tbb_init = new tbb::task_scheduler_init();
}

void CellChain::initD()
{
    double D_san = 0.05;
    double D_atr  = 0.2 ;
    double D_avn = 0.05;
    double D_pur = 0.8;
    double D_vent = 0.2;
    for (int i=0; i<size; i++){
        if (i < sanLen){
            double value = 0.075+0.005*(double)i/((double)sanLen-1);
            ((luo_rudy_I_model_1991 *)cells[i])->setGsi(0.0915);
            D[i] = D_san;
            ((luo_rudy_I_model_1991 *)cells[i])->setGK1(value);
        }
        else if (i < sanLen+atriaLen) {
            D[i] = D_atr;
            ((luo_rudy_I_model_1991 *)cells[i])->setGK1(0.2);
        }
        else if (i < sanLen+atriaLen+avLen){
            ((luo_rudy_I_model_1991 *)cells[i])->setGsi(0.09304);
            double value = 0.08+0.005*((double)(i % (sanLen+atriaLen)))/((double)avLen-1);
            ((luo_rudy_I_model_1991 *)cells[i])->setGK1(value);
            D[i] = D_avn;
        }
        else if (i < sanLen+atriaLen+avLen+purkLen) {
            D[i] = D_pur;
            ((luo_rudy_I_model_1991 *)cells[i])->setGk(0.2);//0.282
        }
        else                                        D[i] = D_vent;

        //D[i] = 0.1;
    }
    int delta = 10;
    for (int i = 0; i<delta; i++){
        D[sanLen-delta/2+i]                       = D_san+(D_atr-D_san)*(double)i/(double)delta;
        D[sanLen+atriaLen-delta/2+i]              = D_atr+(D_avn-D_atr)*(double)i/(double)delta;
        D[sanLen+atriaLen+avLen-delta/2+i]        = D_avn+(D_pur-D_avn)*(double)i/(double)delta;
        D[sanLen+atriaLen+avLen+purkLen-delta/2+i]= D_pur+(D_vent-D_pur)*(double)i/(double)delta;
    }

    FILE *ofs = fopen("couplings.m","w");
    fprintf(ofs,"clear all;\nx=[\n");
    for (int i=0; i<size; i++){
        fprintf(ofs,"%g\n",D[i]);
    }
    fprintf(ofs,"];\nplot(x,'bo-');\n");
    fclose(ofs);
}

/*
  LR_regular: prop_start=0.1 -slow; 1 - fast;
  */

Cell* CellChain::getCell(int i){
    return this->cells[i];
}

void CellChain::SolveDt(double dt, double time){
#if 0
    for (int i=0; i<this->size; i++){
        Solver->IntegrateOverDt(cells[i],cells[i]->dt,time);
    }
    for (int i=0; i<this->size; i++){
        Coupling(i);
    }
#endif
    tbb::parallel_for(tbb::blocked_range<size_t>(0,this->size),
                      ParallelSolver(this,dt,time,0));
    tbb::parallel_for(tbb::blocked_range<size_t>(0,this->size),
                      ParallelSolver(this,dt,time,1));
}

void CellChain::Coupling(int i){
    if (i == 0){
        cells[0]->incrementV((cells[1]->getV()-cells[0]->getV())*cells[0]->dt*this->D[0]);
    } else if (i == size-1){
        cells[size-1]->incrementV((cells[size-2]->getV()-cells[size-1]->getV())*cells[size-1]->dt*this->D[size-1]);
    } else {
        cells[i]->incrementV((cells[i-1]->getV()+cells[i+1]->getV()-2*cells[i]->getV())*cells[i]->dt*this->D[i]);
    }
}


void CellChain::Solve(double MaxTime, double skipTime){
    double dt = cells[0]->dt;
    int saveStep = (int)5/dt;
    int MT = (int)(MaxTime/dt);
    int skip = (int)(skipTime/dt);

    //int avSteps = 5;
    FILE *ofs = fopen("rst.m","w");
    fprintf(ofs,"clear all;\nx=[\n");
    for (int i=0; i<MT; i++){
        SolveDt(dt,i*dt);
        if (i/saveStep*saveStep == i && i>skip){
            printf("%d steps out of %d done\n",i,MT);
            for (int j=0; j<size; j++){
                fprintf(ofs,"%g ",cells[j]->getV());
            }
            fprintf(ofs,"\n");
            fflush(stdout);
        }
    }
    fprintf(ofs,"];\nimagesc(x');\ncolorbar;\n");
    fclose(ofs);
}

static bool inVector(std::vector<int> &vec, int &value){
    for (int i=0; i<vec.size(); i++){
        if (value == vec[i]) return true;
    }
    return false;
}

void CellChain::setStimAmplitudeColl(double A, std::vector<int> ids){
    for (int i=0; i<size; i++){
        if (!ids.size() || inVector(ids,i))
            *(cells[i]->stimulus.stim_amplitude) = A;
    }
}

void CellChain::setStimStartColl(double time){
    for (int i=0; i<size; i++){
        *(cells[i]->stimulus.stim_start) = time;
    }
}

