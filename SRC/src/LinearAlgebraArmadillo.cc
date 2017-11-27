#ifndef _LINEARALGEBRAARMADILLO_CC_
#define _LINEARALGEBRAARMADILLO_CC_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include "LinearAlgebraArmadillo.h"

LinearAlgebraArmadillo::LinearAlgebraArmadillo() : LinearAlgebraBase() {

}

void LinearAlgebraArmadillo::GenerateVSUMatrixAndCoefficientsVector(\
        const vector<vector<double> > &m,const vector<double> &slp, const int J,\
        vector<vector<double> > &vsu,vector<double> &c) {
    arma::mat M;
    M=CopyFromStdMat(m);
    arma::mat U,V;
    arma::vec Sv;
    if (!(arma::svd(U, Sv, V, M))) {
        throw string("Error in svd decomposition (armadillo)!");
    }
    int nr=M.n_rows;
    nr-=(J);
    M.clear();
    arma::mat res(J,nr);
    res.zeros();
    arma::mat R=arma::join_rows(arma::pinv(arma::diagmat(Sv)),res);
    res.clear();
    arma::mat Ut=U.t();
    U.clear();
    arma::mat RUt=R*Ut;
    R.clear();
    arma::mat VSU=(V*RUt);
    V.clear();
    RUt.clear();
    Sv.clear();
    Sv=CopyFromStdVec(slp);
    vsu=CopyToStdMat(VSU);
    c=CopyToStdVec(Sv);
    VSU.clear();
    Sv.clear();
}

vector<double> LinearAlgebraArmadillo::MatrixDotVector(const vector<vector<double> >&m,\
        const vector<double> &v) {
    arma::mat M;
    M=CopyFromStdMat(m);
    arma::vec V,resarma;
    V=CopyFromStdVec(v);
    resarma=M*V;
    M.clear();
    V.clear();
    vector<double> resvec=CopyToStdVec(resarma);
    resarma.clear();
    return resvec;
}

arma::vec LinearAlgebraArmadillo::CopyFromStdVec(const vector<double> &vin) {
    int nn=int(vin.size());
    arma::vec res(nn);
    for ( int i=0 ; i<nn ; ++i ) {
        res[i]=vin[i];
    }
    return res;
}

arma::mat LinearAlgebraArmadillo::CopyFromStdMat(const vector<vector<double> > &min) {
    int nr=int(min.size());
    if ( nr<1 ) {
        cout << "Matrix has zero rows!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
        arma::mat res(1,1);
        return res;
    }
    int nc=int(min[0].size());
    if ( nc<1 ) {
        cout << "Matrix has zero cols!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
        arma::mat res(1,1);
        return res;
    }
    arma::mat res(nr,nc);
    for ( int i=0 ; i<nr ; ++i ) {
        for ( int j=0 ; j<nc ; ++j ) { res(i,j)=min[i][j]; }
    }
    return res;
}

vector<double> LinearAlgebraArmadillo::CopyToStdVec(const arma::vec &vin) {
    size_t nn=size_t(vin.n_elem);
    vector<double> res;
    if ( nn<1 ) {
        cout << "Vector has zero elements!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
        return res;
    }
    res.resize(nn);
    for ( int i=0 ; i<nn ; ++i ) {
        res[i]=vin[i];
    }
    return res;
}

vector<vector<double> > LinearAlgebraArmadillo::CopyToStdMat(const arma::mat &min) {
    vector<vector<double> > res;
    size_t nr=size_t(min.n_rows);
    if ( nr<1 ) {
        cout << "Matrix has zero rows!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
        return res;
    }
    size_t nc=size_t(min.n_cols);
    if ( nc<1 ) {
        cout << "Matrix has zero cols!" << endl;
        cout << __FILE__ << ", line: " << __LINE__ << endl;
        return res;
    }
    res.resize(nr);
    for ( size_t i=0 ; i<nr ; ++i ) { res[i].resize(nc); }
    for ( size_t i=0 ; i<nr ; ++i ) {
        for ( int j=0 ; j<nc ; ++j ) { res[i][j]=min(int(i),int(j)); }
    }
    return res;
}

#endif  /* _LINEARALGEBRAARMADILLO_CC_ */

