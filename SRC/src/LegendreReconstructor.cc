#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <ctime>
#include "LegendreReconstructor.h"
#include "LinearAlgebraArmadillo.h"

// CONSTRUCTOR =================================================================

LegendreReconstructor::LegendreReconstructor(int pl_order, DATA* _data) : WavefrontReconstructor(_data) {
    linalg=new LinearAlgebraArmadillo();
    FindWaveFrontCoefficientsFromCoordinatesAndSlopes(pl_order);
    //return;
}

LegendreReconstructor::~LegendreReconstructor() {
    delete linalg;
    linalg=NULL;
}

// =============================================================================

void LegendreReconstructor::FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int pl_order) {
    N=pl_order;
    J=N*N;

    SetSlopes(data->dx,data->dy);

    // TODO: reimplement timing, or change by runtime profiling
    clock_t c_init,c_end;
    if ( !haveMat ) { // TODO: change syntax getIfhaveMat??
        c_init=clock();
        GenerateMatrixM(data->x,data->y);
        GenerateRHSOfLeastSquareEquation();
        GenerateReconstructionMatrix(data->x,data->y);
        c_end=clock();
        cpuTimeMatGen=double(c_end-c_init)/double(CLOCKS_PER_SEC);
    }

    c_init=clock();
    ComputeCoefficients();
    c_end=clock();

    cpuTimeRHSLSE=double(c_end-c_init)/double(CLOCKS_PER_SEC);
}


// PRIVATE ---------------------------------------------------------------------

void LegendreReconstructor::GenerateMatrixM(const vector<double> &x,const vector<double> &y) {
    size_t nn=x.size();
    int n1,n2;
    double ooN=1.0e0/double(N);
    _M.resize(2*nn);
    for ( size_t i=0 ; i<(2*nn) ; ++i ) {
        _M[i].resize(J);
    }
    for ( int j=0 ; j<J ; ++j ) {
        n1=floor(double(j)*ooN);
        n2=j%N;
        for ( int i=0 ; i<nn ; ++i ) {
            _M[i][j]=(LegendrePolynomials::FirstDerivativeNormalPn(n1,x[i]))*(LegendrePolynomials::NormalPn(n2,y[i]));
            _M[i+nn][j]=(LegendrePolynomials::NormalPn(n1,x[i]))*(LegendrePolynomials::FirstDerivativeNormalPn(n2,y[i]));
        }
    }
}

// -----------------------------------------------------------------------------
void LegendreReconstructor::GenerateReconstructionMatrix( const vector<double> &x, const vector<double> &y) {
    if ( !haveMat ) {
        cout << "Error: First generate the matrix M, and compute\n"
            "coefficients!" << endl;
        return;
    }
    int nn=x.size();
    int jj=J;
    _phases.resize(nn);
    _R.resize(nn);
    for ( int i=0 ; i<nn ; ++i ) {
        _R[i].resize(J);
        for ( int j=0 ; j<jj ; ++j ) { _R[i][j]=0.0e0; }
    }
    double ooN=1.0e0/double(N);
    int n1,n2;
    for ( int j=0 ; j<jj ; ++j ) {
        for ( int i=0 ; i<nn ; ++i ) {
            n1=floor(j*ooN);
            n2=j%N;
            _R[i][j]=((LegendrePolynomials::NormalPn(n1,x[i]))*(LegendrePolynomials::NormalPn(n2,y[i])));
        }
    }
}
