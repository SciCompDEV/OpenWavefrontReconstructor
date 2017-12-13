#ifndef _AUXILIARYTEMPORARYFUNCTIONS_H_
#define _AUXILIARYTEMPORARYFUNCTIONS_H_

#include <iomanip> // setprecision
#include <string>
using std::string;
#include <vector>
using std::vector;
class AuxiliaryTemporaryFunctions {
public:
    static double GetCoefficientOfDetermination(vector<double> &known,vector<double> &estim);
    static double GetCoefficientOfDetermination(vector<double> &known,vector<double> &estim, double &sumtot,double &sumres);
    static double GetNormalizedRMS(vector<double> &known,vector<double> &estim);
    static void SaveData3D(const vector<double> &xx,const vector<double> &yy,const vector<double> &zz,const string fileName=string("temp3d.tsv"));
    static void MakePlot3DGnuplot(string gnuplotZLabel,string gnuplotTitle, string fdat=string("temp3d.tsv"),string fpdf=string("temp3d.pdf"));
    static void MakePlot2DGnuplot(string gnuplotYLabel,string gnuplotTitle, string fdat=string("temp2d.dat"),string fpdf=string("temp2d.pdf"));
    static void SaveComparedData3D(vector<double> &xx,vector<double> &yy,vector<double> &zor,vector<double> &zrec, string fileName=string("temp3dcomp.tsv"));
    static void SaveDiffData2D(const vector<double> &zor,const vector<double> &zrec, string fileName=string("temp2d.dat"));
    // Final gather
    static void SaveData3D_gathered(vector<double>&, vector<double>&, vector<double>&, vector<double>&, string, string, string, string);
    static void print_simu_param_gather(vector<double>&,vector<double>&, std::string, int, int);
};

#endif
