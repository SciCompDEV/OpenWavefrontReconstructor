#ifndef _LINEARALGEBRABASE_H_
#define _LINEARALGEBRABASE_H_

#include <vector>
using std::vector;
class LinearAlgebraBase {
public:
    LinearAlgebraBase() {}
    virtual ~LinearAlgebraBase() {};
    /** As the name suggests, it generates the matrix VSU (see main paper),
     * and the coefficients vector (see vector \f$\boldsymbol A\f$ in the
     * main paper. It requires the matrix m, which contains the derivatives
     * of the polynomials at each point of the grid (see the matrix
     * \f$\boldsymbol M\f$ of the main paper). This is the function that MUST
     * be implemented in all derived classes, and it is also the function
     * used in WavefrontReconstructor class. */

    virtual void GenerateVSUMatrixAndCoefficientsVector(const vector<vector<double> > &m,\
           const vector<double> &slp, const int J,\
            vector<vector<double> > &vsu,vector<double> &c) = 0;
    /** Returns a vector which is the matrix.dot.vector product. This function
     * is used in the WavefrontReconstructor class.  */
    virtual vector<double> MatrixDotVector(const vector<vector<double> > &m,\
            const vector<double> &v) = 0;
protected:
    /* In derived classes, there should be copy functions, from std::vector<double>
     * and std::vector<std::vector<double> > to the specific library
     * vector and matrix containers. For instance, if the library to be
     * used is armadillo, there should be the functions
     * void CopyFromStdVec(const vector<double> &vin,arma::vec &vout);
     * void CopyToStdVec(const arma::vec &vin,vector<double> &vout);
     * These functions cannot be declared virtual, because of the specific
     * types (e.g. arma::vec), but there should be present in the derived
     * classes. */
};


#endif  /* _LINEARALGEBRABASE_H_ */

