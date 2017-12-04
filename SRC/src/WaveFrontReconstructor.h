#ifndef __WavefrontReconstructor_H_
#define __WavefrontReconstructor_H_
#include <vector>
using std::vector;
#include <armadillo>
#include "DATA.h"
#include "LinearAlgebraBase.h"
#include "CLogging.h"
#include <memory>

class WavefrontReconstructor {

protected:
    // TODO: document this ... 
    int J,N;
    vector<double> _slopes; /*!< Gradient vector (see vector \f$\boldsymbol G\f$ in the main paper.  */
    vector<double> _phases; /*!< Internal vector where to store the phases, i.e. the wavefront values  */
    vector<vector<double> > _M;/*!< The derivatives matrix, see matrix \f$\boldsymbol M\f$ of the main paper. */
    vector<vector<double> > _R; /*!< The reconstruction matrix (similar to M, but for the wavefront; M is for gradients.)  */
    vector<vector<double> > _VSU; /*!< The VSU matrix (see the main paper), used to recover coefficients of the
                                    wavefront expansion.  */
    bool haveMat;
    // TODO: documet this
    double cpuTimeMatGen,
           cpuTimeRHSLSE,
           cpuTimeSimpleReconstruction,
           cpuTimeGenMatrixM,
           cpuTimeGenMatrixR,
           cpuTimePureSVD;
    unique_ptr<CLogging> WFREC;

public:

    DATA* data;
    LinearAlgebraBase* linalg;
    vector<double> _coeffs;
    // 
    WavefrontReconstructor(DATA* _data); 
    //
    virtual ~WavefrontReconstructor() {
       ; 
    }

    // -------------------------------------------------------------------------
    // For determining the coefficients of the expansion.
    // to be reimplemented within each particular concrete class 
    virtual void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder) = 0;

    // -------------------------------------------------------------------------
    //  Once the function GenerateReconstructionMatrix, and FindWaveFrontCoefficientsFromCoordinatesAndSlopes had been called,
    //  the coefficients would be simply computed with this function.
    virtual void ComputeCoefficients(void);

    // -------------------------------------------------------------------------
    // Makes the reconstruction using the data contained in _coeffs: will compute the wavefront (and will store the results in "z"
    // Assumes that the estimation of the coefficients was performed before calling it 
    // (the estimation is called by FindWaveFrontCoefficientsFromCoordinatesAndSlopes). 
    virtual vector<double> ComputeReconstructedWaveFront(
            const vector<double> &x,
            const vector<double> &y,
            const vector<double> &coeffs);

    // -------------------------------------------------------------------------
    // The algorithm is in essence an integration -> the wavefront is determined up to an arbitrary constant. 
    // This function remedies somewhat this by centering the resultant wavefront according to the average of zz. 
    virtual void CenterWavefrontAlongZ(vector<double> &zz);

    //Copies the user slopes (possibly from a CCD camera) to the
    //internal slopes data array (in this implementation an arma:vec)
    void SetSlopes(const vector<double> &dx,const vector<double> &dy);

    // -------------------------------------------------------------------------
    // TODO: not as the name suggest
    // CPU time used for computing  all the matrices involved in the process (in seconds).  */
    double GetCPUTimeMatrixGeneration(void) {
        return cpuTimeMatGen;
    }

    // -------------------------------------------------------------------------
    // Return the CPU time employed to compute the coefficients 
    double GetCPUTimeCoefficientEstimation(void) {
        return cpuTimeRHSLSE;
    }
    // Returns the CPU time spent on performing a simple reconstruction,
    // i.e. the product _R*_coeffs
    double GetCPUTimeSimpleReconstruction(void) {
        return cpuTimeSimpleReconstruction;
    }
    //Returns the CPU time spent on assembling the matrix M (see Eq.(8) of
    //the main paper.
    double GetCPUTimeGenerationOfMatrixM() { return cpuTimeGenMatrixM; }
    //Returns the CPU time spent on assembling the matrix R (see Eq.(11) of
    //the main paper.
    double GetCPUTimeGenerationOfMatrixR() { return cpuTimeGenMatrixR; }
    //Returns the CPU time spent on performing exclusively the SVD operation.
    double GetCPUTimePureSVD() { return cpuTimePureSVD; }

    // -------------------------------------------------------------------------
    // Returns the actual number of terms used in the expansion.
    int NumberOfPolynomialTerms(void) { return J; }
    // Returns the polynomial order (max n of Ynm).
    int PolynomialOrder(void) { return N; }

    // -------------------------------------------------------------------------
    void PrintCoefficients(void);
    virtual string PolynomialType()=0;


   virtual  void print_CPU_gathered_info();

protected:

    //
    // To generate the matrix that will be solved through the SVD algorithm = the matrix with the derivatives of the polynomials at the pairs (x,y). 
    virtual void GenerateMatrixM(const vector<double> &x,const vector<double> &y) = 0;

    // To generate the matrix used for the reconstruction, after th Least-Squares method was applied. It contains the values of the polynomials at the points (x,y). 
    virtual void GenerateReconstructionMatrix(const vector<double> &x, const vector<double> &y) = 0;

    // For solving the reconstruction problem. Usually there would be no need to re-implement this function.  
    virtual void GenerateRHSOfLeastSquareEquation(void);

    //  
    //void CopyWholeArray(const vector<double> &vin,arma::vec &vout, size_t posStart,size_t nElem);

};


#endif  /* __WavefrontReconstructor_H_ */
