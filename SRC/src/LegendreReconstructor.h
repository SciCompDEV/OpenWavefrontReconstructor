#ifndef LEGENDRERECONSTRUCTOR_H_
#define LEGENDRERECONSTRUCTOR_H_
#include "WaveFrontReconstructor.h"
#include "LegendrePolynomials.h"
#include "DATA.h"
class LegendreReconstructor : public WavefrontReconstructor {
    public:
        string PolynomialType() {return string("Legendre");}
    protected:
        LegendrePolynomials lp;

    public:
        LegendreReconstructor(int pl_order, DATA* _data); 
        ~LegendreReconstructor();
        // See wavefrontreconstructor.h for more information
        void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder);
    protected:

        // See wavefrontreconstructor.h for more information 
        void GenerateMatrixM(const vector<double> &x,const vector<double> &y);

        // See wavefrontreconstructor.h for more information  
        void GenerateReconstructionMatrix(const vector<double> &x, const vector<double> &y);

};
#endif  /* LEGENDRERECONSTRUCTOR_H_ */
