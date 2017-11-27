#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <ctime>
#include "HalfCircularHarmonicsReconstructor.h"
#include "LinearAlgebraArmadillo.h"

// CONSTRUCTOR =================================================================

HalfCircularHarmonicsReconstructor::HalfCircularHarmonicsReconstructor(int pl_order, DATA* _data) : WavefrontReconstructor(_data) {
    linalg=new LinearAlgebraArmadillo();
    FindWaveFrontCoefficientsFromCoordinatesAndSlopes(pl_order);
    //return;
}

HalfCircularHarmonicsReconstructor::~HalfCircularHarmonicsReconstructor() {
    delete linalg;
    linalg=NULL;
}

// =============================================================================

void HalfCircularHarmonicsReconstructor::FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int pl_order) {
    N=pl_order;
    J=hch.ComputeNumberOfPolynomialsFromPolynomialOrder(N);
    SetSlopes(data->dx,data->dy);

    // TODO: reimplement timing, or change by runtime profiling
    clock_t c_init,c_end,c_tmp1;
    if ( !haveMat ) { // TODO: change syntax getIfhaveMat??
        c_init=clock();

        c_tmp1=clock();
        GenerateMatrixM(data->x,data->y);
        cpuTimeGenMatrixM=double(clock()-c_tmp1)/double(CLOCKS_PER_SEC);

        c_tmp1=clock();
        GenerateRHSOfLeastSquareEquation();
        cpuTimePureSVD=double(clock()-c_tmp1)/double(CLOCKS_PER_SEC);

        c_tmp1=clock();
        GenerateReconstructionMatrix(data->x,data->y);
        cpuTimeGenMatrixR=double(clock()-c_tmp1)/double(CLOCKS_PER_SEC);
        c_end=clock();

        cpuTimeMatGen=double(c_end-c_init)/double(CLOCKS_PER_SEC);
    }

    c_init=clock();
    ComputeCoefficients();
    c_end=clock();

    cpuTimeRHSLSE=double(c_end-c_init)/double(CLOCKS_PER_SEC);
}


// PRIVATE ---------------------------------------------------------------------

void HalfCircularHarmonicsReconstructor::GenerateMatrixM(const vector<double> &x,const vector<double> &y) {
    size_t nn=x.size();
    _M.resize(2*nn);
    for ( size_t i=0 ; i<(2*nn) ; ++i ) {
        _M[i].resize(J);
    }
    for ( int j=1 ; j<=J ; ++j ) {
        for ( size_t i=0 ; i<nn ; ++i ) {
            _M[i][j-1]=hch.DCircularHarmonicDxiFromCartesian(j,x[i],y[i]);
            _M[i+nn][j-1]=hch.DCircularHarmonicDetaFromCartesian(j,x[i],y[i]);
        }
    }
}

// -----------------------------------------------------------------------------
void HalfCircularHarmonicsReconstructor::GenerateReconstructionMatrix( const vector<double> &x, const vector<double> &y) {
    if ( !haveMat ) {
        cout << "Error: First generate the matrix M, and compute\n"
             "coefficients!" << endl;
        return;
    }
    int nn=x.size();
    int jj=J;
    _phases.resize(nn);
    _R.resize(nn);
    for ( int i=0 ; i<nn ; ++i ) { _R[i].resize(J); }
    for ( int j=0 ; j<jj ; ++j ) {
        for ( int i=0 ; i<nn ; ++i ) {
            _R[i][j]=hch.HalfCircularHarmonicFromCartesian((j+1),x[i],y[i]);
        }
    }
}
