#ifndef _LINEARALGEBRAARMADILLO_H_
#define _LINEARALGEBRAARMADILLO_H_

#include "LinearAlgebraBase.h"
#include <armadillo>

class LinearAlgebraArmadillo : public LinearAlgebraBase {
public:
    LinearAlgebraArmadillo();
    void GenerateVSUMatrixAndCoefficientsVector(
            const vector<vector<double> > &m,
            const vector<double> &slp, 
            const int J,
            vector<vector<double> > &vsu,
            vector<double> &c);
    vector<double> MatrixDotVector(
            const vector<vector<double> > &m,
            const vector<double> &v);
    static arma::vec CopyFromStdVec(const vector<double> &vin);
    static arma::mat CopyFromStdMat(const vector<vector<double> > &min);
    static vector<double> CopyToStdVec(const arma::vec &vin);
    static vector<vector<double> > CopyToStdMat(const arma::mat &min);
protected:
};


#endif  /* _LINEARALGEBRAARMADILLO_H_ */

