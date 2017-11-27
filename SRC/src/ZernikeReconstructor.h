#ifndef ZERNIKEWFRECONSTRUCTOR_H_
#define ZERNIKEWFRECONSTRUCTOR_H_
#include "WaveFrontReconstructor.h"
#include "ZernikePolynomials.h"
#include "DATA.h"
class ZernikeReconstructor : public WavefrontReconstructor {
    public:
        string PolynomialType() {return string("CircularZernike");}
    protected:
        ZernikePolynomials zp;

    public:
        ZernikeReconstructor(int pl_order, DATA* _data); 
        ~ZernikeReconstructor();
        // See wavefrontreconstructor.h for more information
        void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder);
    protected:

        // See wavefrontreconstructor.h for more information 
        void GenerateMatrixM(const vector<double> &x,const vector<double> &y);

        // See wavefrontreconstructor.h for more information  
        void GenerateReconstructionMatrix(const vector<double> &x, const vector<double> &y);

};
#endif  /* ZERNIKEWFRECONSTRUCTOR_H_ */
