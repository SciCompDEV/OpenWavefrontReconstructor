#ifndef _HALFCIRCULARHARMONICS_H_
#define _HALFCIRCULARHARMONICS_H_

class HalfCircularHarmonics {
public:
    /** Returns the Associated Legendre Polynomial of order n,m, at x.  */
    static double NormalPnm(int n,int m,double x);
    /** Returns the first derivative of Pnm, at x.  */
    static double FirstDerivativeNormalPnm(int n,int m,double x);
    /** Returns the first derivative of Pn,, atx, but using previously 
        calculated Pnm (which must be passed as arguments). */
    static double FirstDerivativeNormalPnmFromPnmAndPnplus1m(int n,int m,\
            double x,double pnm,double pnp1m);
    /** Returns the first derivative of Pn,, atx, but using previously 
        calculated Pnm (which must be passed as arguments). */
    static double FirstDerivativeNormalPnmFromPnmAndPnminus1m(int n,int m,\
            double x,double pnm,double pnm1m);
    /** Returns the order of Pnm (i.e., n and m) given a one
        index ordered polynomial.  */
    static inline void ComputeNMFromJ(int j,int &n,int &m){
        int lm=j-1, ln=0;
        while ( lm>ln ) {++ln; lm-=(2*ln);}
        n=ln; m=lm;
    }
    /** Given the principal order n, the function returns
        the number of polynomials that forms a closed shell.  */
    static int ComputeNumberOfPolynomialsFromPolynomialOrder(int n);
    static constexpr double cOneOverSqrtPi=0.564189583547756286948079e0; //1/sqrt{pi}
    static constexpr double cOneOverSqrtTwoPi=0.398942280401432677939946e0; // 1/sqrt{2pi}
    static constexpr double cPi=3.14159265358979323846264e0;
    static constexpr double cHalfPi=1.57079632679489661923132e0;
    /** Computes \f$\mu\f$ and \f$\varphi\f$ from \f$\xi\f$ and \f$\eta\f$ (see main paper).
        This is, the composite mapping \f$\chi_3\circ\chi_2\circ\chi_1\f$. The computed \f$\mu\f$ and \f$\phi\f$
        are saved in the parameters mu and phi. */
    static void ComputeThetaMuAndPhiFromXiEta(const double xi,const double eta,double &theta,double &mu,double &phi);
    static double DCircularHarmonicDxiFromCartesian(const int n,const int m,const double xi,const double eta);
    static inline double DCircularHarmonicDxiFromCartesian(const int j,const double xi,const double eta) {
        int n,m; ComputeNMFromJ(j,n,m); return DCircularHarmonicDxiFromCartesian(n,m,xi,eta);
    }
    static double DCircularHarmonicDxi(const int n,const int m,const double theta,
            const double mu,const double phi);
    static double DCircularHarmonicDetaFromCartesian(const int n,const int m,const double xi,const double eta);
    static inline double DCircularHarmonicDetaFromCartesian(const int j,const double xi,const double eta) {
        int n,m; ComputeNMFromJ(j,n,m); return DCircularHarmonicDetaFromCartesian(n,m,xi,eta);
    }
    static double DCircularHarmonicDeta(const int n,const int m,const double theta,
            const double mu,const double phi);
    /* ************************************************************************** */
    static double CircularHarmonic(const int n,const int m,const double mu,const double phi);
    static double HalfCircularHarmonicFromCartesian(int n,int m,double x,double y);
    static inline double HalfCircularHarmonicFromCartesian(int j,double x,double y){
        int n,m; ComputeNMFromJ(j,n,m); return HalfCircularHarmonicFromCartesian(n,m,x,y);
    }
};


#endif  /* _HALFCIRCULARHARMONICS_H_ */

