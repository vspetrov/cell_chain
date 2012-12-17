#ifndef GNUPLOTPIPER_H
#define GNUPLOTPIPER_H
#include <stdio.h>
#include <vector>
using namespace std;

class GnuplotPiper{
private:
    double dt;
    unsigned short skipFactor;
    int startPosition;
    FILE *gnuplotPipe;
    vector<double> time;
    vector<double> value;
    unsigned long int replotSkip;
public:
    GnuplotPiper(unsigned short ReplotNumber, double MaxTime = 1000, double Dt = 0.005, unsigned short SkipFactor = 1000){
        gnuplotPipe = popen("gnuplot --persist","w");
        fprintf(gnuplotPipe,"set yrange[-100:50]\n");
        this->startPosition = 0;
        this->replotSkip = ((unsigned long int)(MaxTime/Dt))/ReplotNumber;
        this->dt = Dt;
        this->skipFactor = SkipFactor;
    }
    ~GnuplotPiper(){
        fprintf(gnuplotPipe,"exit \n");
        pclose(gnuplotPipe);
    }
    void pipe1D(long unsigned int &intTime, double value){
        if (intTime/skipFactor*skipFactor == intTime){
            this->time.push_back(intTime*dt);
            this->value.push_back(value);
        }
        if (intTime/replotSkip*replotSkip == intTime)
            plot();
    }
private:
    void plot(){
        reinitGnuplot();
        for (unsigned int i=0; i<time.size(); i++){
            fprintf(gnuplotPipe,"%g\t%g\n",time.at(i),value.at(i));
        }
        fprintf(gnuplotPipe,"e\n");
        fflush(gnuplotPipe);
    }
    void reinitGnuplot(){
        fprintf(gnuplotPipe,"plot '-' w l\n");
    }
};

#endif // GNUPLOTPIPER_H
