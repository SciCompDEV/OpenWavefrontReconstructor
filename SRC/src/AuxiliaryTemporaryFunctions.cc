#ifndef _AUXILIARYTEMPORARYFUNCTIONS_CC_
#define _AUXILIARYTEMPORARYFUNCTIONS_CC_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <cmath>
#include "AuxiliaryTemporaryFunctions.h"

double AuxiliaryTemporaryFunctions::GetCoefficientOfDetermination(vector<double> &known,
        vector<double> &estim) {
    double st,sr;
    return GetCoefficientOfDetermination(known,estim,st,sr);
}

double AuxiliaryTemporaryFunctions::GetCoefficientOfDetermination(vector<double> &known,\
      vector<double> &estim,double &st,double &sr) {
   int nn=known.size();
   if ( nn!=int(estim.size()) ) {
      cerr << "Warning: vectors have not the same size" << endl;
      cerr << "file: " << __FILE__ << ", line: " << __LINE__ << endl;
      if ( nn<int(estim.size()) ) {
         nn=int(estim.size());
      }
   }
   double tot=0.0e0;
   double res=0.0e0;
   double av=0.0e0;
   for ( int i=0 ; i<nn ; ++i ) { av+=known[i]; }
   av/=double(nn);
   double tmp;
   for ( int i=0 ; i<nn ; ++i ) {
      tmp=known[i]-av;
      tot+=(tmp*tmp);
      tmp=known[i]-estim[i];
      res+=(tmp*tmp);
   }
   st=tot;
   sr=res;
   return (1.0e0-(res/tot));
}

double AuxiliaryTemporaryFunctions::GetNormalizedRMS(vector<double> &known,vector<double> &estim) {
   int nn=known.size();
   if ( nn!=int(estim.size()) ) {
      cerr << "Warning: vectors have not the same size" << endl;
      cerr << "file: " << __FILE__ << ", line: " << __LINE__ << endl;
      if ( nn<int(estim.size()) ) {
         nn=int(estim.size());
      }
   }
    double diff=0.0e0,knsq=0.0e0;
    for ( int i=0 ; i<nn ; ++i ) {
        diff+=((known[i]-estim[i])*(known[i]-estim[i]));
        knsq+=(known[i]*known[i]);
    }
    return (diff/knsq);
}

void AuxiliaryTemporaryFunctions::SaveData3D(const vector<double> &xx,const vector<double> &yy,const vector<double> &zz,const string fileName) {
   ofstream dfil(fileName.c_str());
   //dfil << "#Number of maximum points per direction (used to build the grid): "
        //<< NUM_OF_POINTS_PER_DIRECTION << endl;
   //dfil << "#Polynomial type: " << BASE_POLYNOMIALS << endl;
   //dfil << "#Polynomial order: " << POLYNOMIAL_ORDER << endl;
   //dfil << "#Absolute number of polynomials: " << absNofPolTerms << endl;
   //dfil << "#CPU time to perform matrix decomposition: " << (1000.0e0*timeMatOps)
        //<< " milliseconds" << endl;
   //dfil << "#CPU time to reconstruct wavefront: " << (1000.0e0*timeRecWF)
        //<< " milliseconds"<< endl;
   //dfil << "#Mock wavefront type: " << TYPE_OF_MOCK_WAVEFRONT << endl;
   //dfil << << endl;
   int nn=xx.size();
   //cout << "nn: " << nn << endl;
   double tmp=xx[0];
   for ( int i=0 ; i<nn ; ++i ) {
      if ( xx[i]!=tmp ) {
         dfil << endl;
         tmp=xx[i];
      }
      dfil << xx[i] << " " << yy[i] << " " << zz[i] << endl;
   }
   dfil.close();
}
void AuxiliaryTemporaryFunctions::SaveComparedData3D(vector<double> &xx,vector<double> &yy, vector<double> &zor,
        vector<double> &zrec,string fileName) {
   ofstream dfil(fileName.c_str());
   //dfil << "#Number of maximum points per direction (used to build the grid): "
        //<< NUM_OF_POINTS_PER_DIRECTION << endl;
   //dfil << "#Polynomial type: " << BASE_POLYNOMIALS << endl;
   //dfil << "#Polynomial order: " << POLYNOMIAL_ORDER << endl;
   //dfil << "#Absolute number of polynomials: " << absNofPolTerms << endl;
   //dfil << "#CPU time to perform matrix decomposition: " << (1000.0e0*timeMatOps)
        //<< " milliseconds" << endl;
   //dfil << "#CPU time to reconstruct wavefront: " << (1000.0e0*timeRecWF)
        //<< " milliseconds"<< endl;
   //dfil << "#Mock wavefront type: " << TYPE_OF_MOCK_WAVEFRONT << endl;
   //dfil << << endl;
   int nn=xx.size();
   //cout << "nn: " << nn << endl;
   double tmp=xx[0];
   for ( int i=0 ; i<nn ; ++i ) {
      if ( xx[i]!=tmp ) {
         dfil << endl;
         tmp=xx[i];
      }
      dfil << xx[i] << " " << yy[i] << " " << (zor[i]-zrec[i]) << endl;
   }
   dfil.close();
}
// TODO: Solano some unusued parameters: yy
void AuxiliaryTemporaryFunctions::SaveDiffData2D(const vector<double> &xx,const vector<double> &yy, const vector<double> &zor,const vector<double> &zrec,string fileName) {
    ofstream dfil(fileName.c_str());
    int nn=xx.size();
    int idx=0;
    double err;
    for ( int i=0 ; i<nn ; ++i ) {
        if ( zor[i]<0.000001 ) {
            err=0.0e0;
        } else {
            err=fabs(1.0e0-(zrec[i]/zor[i]));
            err*=100.0e0;
        }
        dfil << (idx++) << " " << err << endl;
    }
    dfil.close();
}


void AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string gnuplotZLabel,string gnuplotTitle,\
        string fdat,string fpdf) {
    string cmd="gnuplot -e \"set term postscript eps color enhanced fontscale 1.5";
    cmd+="; set xtics -1,0.4,1.0; set xlabel 'x'";
    cmd+="; set ytics -1,0.4,1.0; set ylabel 'y'";
    cmd+="; set xyplane relative 0.1";
    cmd+=string("; set zlabel '")+gnuplotZLabel+string("' rotate by 90");
    cmd+=string(";set title '")+gnuplotTitle+string("' offset character 0,-2");
    cmd+="; set palette rgbformulae 33,13,10";
#if USEGNUPLOTINTERP
    cmd+="; set dgrid3d 75,75 splines";
    cmd+="; set pm3d at s depthorder border lw 0.25 lt rgb '#000000'";
#else
    cmd+="; set pm3d at s depthorder border lw 0.25 lt rgb 'black'";
#endif
    cmd+="; set output '|epstopdf --filter --outfile=";
    cmd+=fpdf;
    cmd+="'; splot '";
    cmd+=fdat;
    cmd+="' u 1:2:3 w pm3d notitle\"";
    system(cmd.c_str());
}

void AuxiliaryTemporaryFunctions::MakePlot2DGnuplot(string gnuplotYLabel,string gnuplotTitle,string fdat,string fpdf) {
    string cmd="gnuplot -e \"set term postscript eps color enhanced fontscale 1.5";
    cmd+=string("; set ylabel '")+gnuplotYLabel+string("'");
    cmd+=string(";set title '")+gnuplotTitle+string("' offset character 0,-2");
    cmd+="; set output '|epstopdf --filter --outfile=";
    cmd+=fpdf;
    cmd+="'; plot '";
    cmd+=fdat;
    cmd+="' u 1:2 w l lw 2 notitle\"";
    system(cmd.c_str());
}

void AuxiliaryTemporaryFunctions::SaveData3D_gathered(vector<double>& x,
        vector<double>& y,
        vector<double>& z,
        vector<double>& zreconstructed,
        string generated,
        string reconstructed,
        string difference,
        string error) {

        SaveData3D(x,y,z, generated);
        SaveData3D(x,y,zreconstructed, reconstructed);
        SaveComparedData3D(x,y,z, zreconstructed, difference);
        SaveDiffData2D(x,y,z,zreconstructed, error);

}

#endif 

