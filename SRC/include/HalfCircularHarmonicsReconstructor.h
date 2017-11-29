#ifndef HALFCIRCULARHARMONICSWFRECONSTRUCTOR_H_
#define HALFCIRCULARHARMONICSWFRECONSTRUCTOR_H_
#include "WaveFrontReconstructor.h"
#include "HalfCircularHarmonics.h"
#include "DATA.h"
class HalfCircularHarmonicsReconstructor : public WavefrontReconstructor {
    public:
        string PolynomialType() {return string("HalfCircularHarmonics");}
    protected:
        //int N; // polinomial number? Re: No (it's the principal order, i.e. the maximum order used for Ynm).
        HalfCircularHarmonics hch;

    public:
        HalfCircularHarmonicsReconstructor(int pl_order, DATA* _data); 
        ~HalfCircularHarmonicsReconstructor();
        // See wavefrontreconstructor.h for more information
        void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder);
    protected:

        // See wavefrontreconstructor.h for more information 
        void GenerateMatrixM(const vector<double> &x,const vector<double> &y);

        // See wavefrontreconstructor.h for more information  
        void GenerateReconstructionMatrix(const vector<double> &x, const vector<double> &y);

};
#endif  /* HALFCIRCULARHARMONICSWFRECONSTRUCTOR_H_ */
