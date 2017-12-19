#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <iomanip>
#include <string>
using std::string;
#include <functional>
#include "WaveFrontReconstructor.h"

// make_unique not yet implemented in C++11
using namespace std;
	template<typename T, typename... Args>
	std::unique_ptr<T> make_unique(Args&&... args) {
		return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}

// CONSTRUCTOR =================================================================

WavefrontReconstructor::WavefrontReconstructor(DATA* _data) {
   WFREC = make_unique<CLogging>( (char*)"Log4cxxConfig.xml","WFREC") ;
   WFREC->INFO((char *)"Building up WavefrontReconstructor  ");

    data=move(_data);
    N=J=0;

    linalg=NULL;

    _slopes.clear();
    _M.clear();
    _R.clear();
    _phases.clear();
    _VSU.clear();
    _coeffs.clear();

    haveMat=false;
    cpuTimeMatGen=0.0e0;
    cpuTimeRHSLSE=0.0e0;
}

// =============================================================================

void WavefrontReconstructor::GenerateRHSOfLeastSquareEquation(void) {
    if ( linalg!=NULL ) {
        linalg->GenerateVSUMatrixAndCoefficientsVector(_M,_slopes,J,_VSU,_coeffs);
    } else {
        cout << "Error: linear algebra library not defined!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
        return;
    }
    haveMat=true;
}

void WavefrontReconstructor::CenterWavefrontAlongZ(vector<double> &zz) {
    int nn=zz.size();
    if ( nn<=0 ) {
        cerr << "Error: non valid size of array!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
    }
    double sum=0.0e0;
    for ( int i=0 ; i<nn ; ++i ) {
        sum+=zz[i];
    }
    sum/=double(nn);
    for ( int i=0 ; i<nn ; ++i ) {
        zz[i]-=sum;
    }
}

void WavefrontReconstructor::SetSlopes(const vector<double> &_dx, const vector<double> &_dy) {
    size_t nn=_dx.size();
    if ( _dx.size() !=_dy.size() ) {
    // TODO: user error class
       cout << "Error: gradient arrays have different sizes!" << endl;
       cout << __FILE__ << ", line: " << __LINE__ << endl;
    }
    if ( _slopes.size()!=(2*nn) ) {
        _slopes.resize(2*nn);
    }
    // TODO: one single array why?
    for ( size_t i=0 ; i<nn ; ++i ) { _slopes[i]=_dx[i]; }
    for ( size_t i=0 ; i<nn ; ++i ) { _slopes[nn+i]=_dy[i]; }
}

void WavefrontReconstructor::ComputeCoefficients(void) {
   if ( linalg!=NULL ) {
       _coeffs=linalg->MatrixDotVector(_VSU,_slopes); // TODO Alejandro the optimization
   } else {
    cout << "No linear algebra library has been chosen!" << endl;
    cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
}

// TODO: short description ...
vector<double> WavefrontReconstructor::ComputeReconstructedWaveFront() {
    vector<double> z; // TODO: this as a class member, to void continuous creation
    // TODO: reimplement timing, or change by runtime profiling
    clock_t c_init,c_end;
    c_init=clock();
    size_t nn=data->x.size();
    if ( (nn)!=(data->y.size()) ) {
        cout << "Error: arrays x and y are of different sizes!!" << endl
             << "File: " << __FILE__ << ", line: " << __LINE__ << endl;
    }
    if ( z.size()!=nn ) {
        z.resize(nn);
    }
    if ( linalg!=NULL ) {
        _phases=linalg->MatrixDotVector(_R,_coeffs);
    } else {
        cout << "Error: no linear algebra has been chosen!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
    }
    if ( nn!=(_phases.size()) ) {
        cout << "Error: the size of the output array is not\n"
             "equal to the internal size of phases!\n"
             "File: " << __FILE__ << ", line: " << __LINE__  << endl;
    }

   for ( size_t i=0 ; i<nn ; ++i ) { z[i]=_phases[i]; }
   c_end=clock();
   cpuTimeSimpleReconstruction=double(c_end-c_init)/double(CLOCKS_PER_SEC);
   return z;
}


// UTILS =======================================================================

void WavefrontReconstructor::PrintCoefficients(void) {
    cout << std::scientific << std::setprecision(12) << endl;
    size_t nn=_coeffs.size();
    for ( size_t i=0 ; i<nn ; ++i ) {
        cout << _coeffs[i] << endl;
    }
}

void WavefrontReconstructor::print_CPU_gathered_info(){ 

        // TODO: move all this stuff to sringstream and pass to the logger

        cout << scientific << setprecision(10) << endl;
        cout << "NumberOfNodes: " << data->dx.size() << endl;
        cout << "CPUTimeSVD: " << GetCPUTimeMatrixGeneration() << endl;
        cout << "CPUTimePureSVD: " << GetCPUTimePureSVD() << endl;
        cout << "CPUTimeGenMatrixM: " << GetCPUTimeGenerationOfMatrixM() << endl;
        cout << "CPUTimeGenMatrixR: " << GetCPUTimeGenerationOfMatrixR() << endl;
        cout << "SingleCPUTimeVSUGProduct: " << GetCPUTimeCoefficientEstimation() << endl;
        cout << "SingleCPUTimeRAProduct: " << GetCPUTimeSimpleReconstruction() << endl;

}

// =============================================================================
