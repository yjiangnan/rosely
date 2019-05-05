//
//  integral.cpp
//  
//
//  Created by Jiang-Nan Yang on 9/10/18.
//
//
#include <math.h>
#include <limits>

#include "integral.hpp"

double fVls2(double V, double s2, double n, double ElogV, double VlogV0, double Vmax){
    double logV = log(V), d = logV-ElogV, logVmax = log(Vmax), dVmax = logVmax - ElogV;
    return exp(-(logV-logVmax) * (n+1.)/2.  - (d*d - dVmax*dVmax) / 2. / VlogV0  -  (n-1.)*s2/2.*(1/V-1/Vmax) );
}

double EVfV  (double V, double s2, double n, double ElogV, double VlogV0, double Vmax, double M){
    return V        * fVls2(V, s2, n, ElogV, VlogV0, Vmax) / M;
}
double EsgmfV(double V, double s2, double n, double ElogV, double VlogV0, double Vmax, double M){
    return sqrt(V)  * fVls2(V, s2, n, ElogV, VlogV0, Vmax) / M;
}
double VsgmfV(double V, double s2, double n, double ElogV, double VlogV0, double Vmax, double M, double Esgm){
    double d = sqrt(V) - Esgm;
    return d * d    * fVls2(V, s2, n, ElogV, VlogV0, Vmax) / M;
}



double fs2lV (double s2,  double n, double V){
    return exp(log((n-1)/2/V)*((n-1)/2) - lgamma((n-1)/2) + log(s2)*((n-3)/2) -(n-1)*s2/(2*V));
}
double fs2alV(double s2a, double n, double V, double a){
    return fs2lV(pow(s2a, 1./a), n, V) / (a * pow(s2a, (a-1.)/a));
}

double  fVls2a(double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax){
    if (n0==0){
        double d = pow(V, a) - EVa, dmax = pow(Vmax, a) - EVa;
        return exp( log(V/Vmax)*(a-n*0.5-1) - (d*d - dmax*dmax)/2./Vga - (n-1.)*s2/2.*(1/V-1/Vmax) );
    }
    else{
        return exp( (log(V/Vmax) * (n0-n-2) - (n-1)*s2*(1./V-1./Vmax) - (n0-1)*(V-Vmax)/Vga) / 2 );
    }
}

double EVfVa  (double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax, double M){
    return V        * fVls2a(V, s2, n, EVa, Vga, a, n0, Vmax) / M;
}
double EsgmfVa(double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax, double M){
    return sqrt(V)  * fVls2a(V, s2, n, EVa, Vga, a, n0, Vmax) / M;
}
double VsgmfVa(double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax, double M, double Esgm){
    double d = sqrt(V) - Esgm;
    return d * d    * fVls2a(V, s2, n, EVa, Vga, a, n0, Vmax) / M;
}


double ms2fs2a(double s2a, double n, double V, double a){
    return s2a * fs2alV(s2a, n, V, a);
}
double Vsmplfs2a(double s2a, double n, double V, double a, double ms2){
    double d = s2a - ms2;
    return d * d * fs2alV(s2a, n, V, a);
}





#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

double erfinv (double x, double p2=0.) // Calculate z-score from p-value
{
    double x2, r, y;
    int  sign_x;
    
    if (x < -1 || x > 1)
        return NAN;
    
    if (x == 0)
        return 0;
    
    if (x > 0)
        sign_x = 1;
    else {
        sign_x = -1;
        x = -x;
    }
    
    if (x <= 0.7) {
        x2 = x * x;
        r = x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
        r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 +
              erfinv_b1) * x2 + erfinv_b0;
    }
    else {
        double x1 = 1 - x;
        if (p2) x1 = p2;
        y = sqrt (-log (x1 / 2));
        r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
        r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
    }
    
    r = r * sign_x;
    x = x * sign_x;
    
    r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
    r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
    if (r != r) r = -std::numeric_limits<double>::infinity(); //NaN
    
    return r;
}

double z_score(double p){
    return sqrt(2.) * erfinv(2*p - 1, 2*p);
}
