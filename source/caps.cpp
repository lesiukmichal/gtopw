#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <complex>
#include <limits>
#include "omp.h"
#include "../headers/gtopw.h"
#include "../headers/caps.h"
#include "../headers/auxfun1.h"

using namespace std;

extern "C" { 
    #include <quadmath.h>
}

void CAPRadFac( const int lmax, const int nmax,
                const int mpar,
                const double k, const double p,
                double** rad_fac ) {
    /* radial factors for the atomic CAP integrals */
    
    if( abs( k ) > 10.0 * numeric_limits<double>::epsilon() ) {
        cout << " This version of CAP code works only for pure Gaussians!" << endl;
        cout << " Emergency halt." << endl;
        exit( EXIT_FAILURE );
    };
    
    double ep;
    int mt;
    ep = exp( -p );
    mt = mpar + nmax;
    
    FnInts <double> ints( mt, p + p, p );
    FnIntsCalc( ints );
    
    for(int n=0; n<=nmax; n++) {
        for(int i=0; i<=n; i++)
            rad_fac[0][n] += binom[i][n] * ints.fn[mpar+i];
        rad_fac[0][n] = rad_fac[0][n] * ep;
    };
    /* */
};

void FnIntsCalc( FnInts <double> & v ) {
    /* calculation of the basic integrals F_n(\alpha,\beta) */
    if( abs( v.b ) < 10.0 * numeric_limits<double>::epsilon() ) {
        v.fn[0] = 1.0 / v.a;
        for(int n=1; n<=v.nmax; n++) 
            v.fn[n] = v.fn[n-1] * n / v.a;
        return ;
    };
    double start,arg,aa,sa,ps,tt,i0;
    int nn,nstart;
    arg   = v.a * v.a / ( 4.0 * v.b );
    sa    = sqrt( arg );
    start = 0.0;
    i0    = 0.5 * sqrt( M_PI / v.b ) * eerfc( sa );
    
    if( v.nmax == 0 ) { v.fn[0] = i0; return ; };
    if( v.nmax > 50 ) {
        cout<<"The value of nmax is too large in FnIntsCalc!"<<endl;
        cout<<"Emergency halt!"<<endl;
        exit( EXIT_FAILURE );
    };
    
    if( arg > 1.0 ) {
        /* for z > 1: the lower n, the faster the convergence of CF */
        nstart = v.nmax;
        nn = ( nstart%2 == 1 ) ? ( nstart + 1 ) / 2 : nstart / 2;
        aa = ( nstart%2 == 1 ) ? 0.0 : 0.5;
        start = hyperg_U_CF1(aa,nn,arg);
        start = start / ( aa + nn );
        start = start * nstart * ( nstart + 1 ) / ( 4.0 * v.b );
    } else {
        /* for z < 1: here we have no choice but to use n = 50 */
        nstart = 50;
        ps     = 1.0;
        start  = 0.0;
        for(int k=0; k<=50; k++) {
            tt     = u0[k] * ps;
            start += tt;
            ps     = ps * sa;
            if( abs( tt / start ) < numeric_limits<double>::epsilon() ) break;
        };
        start = start * nstart * ( nstart + 1 ) / ( 4.0 * v.b );
    };
    double* tn = new double[nstart];
    memset(tn,0,nstart*sizeof(double));
    
    tn[nstart-1] = nstart / v.a - 2.0 * v.b * start / v.a;
    for(int n=nstart-1; n>=1; n--) tn[n-1] = n / ( v.a + 2.0 * v.b * tn[n] );
    
    v.fn[0] = i0;
    for(int n=1; n<=v.nmax; n++) v.fn[n] = tn[n-1] * v.fn[n-1];
    /* */
    delete []tn;
};

double eerfc( const double x ) {
    /*    special code for the scaled error function, exp(x^2) erfc(x)   */
    if ( x < 25.0 ) return exp( x * x ) * erfc( x );
    if ( x < 100.0 ) {
        static double P[] = { 
        +5.641895835477562e-01, -2.820947917738019e-01, +4.231421872959319e-01,
        -1.057854676612069e+00, +3.701638603055076e+00, -1.621315316006683e+01 };
        
        double px2,val;
        px2 = 1.0 / ( x * x );
        val = P[5];
        for(int n=4; n>=0; n--) val = val * px2 + P[n];
        val = val / x;
        return val;
    } else {
        static double R[] = { 
        +5.641895835477563e-01, -2.820947917738781e-01, +4.231421876608172e-01,
        -1.057855469152043e+00, +3.702494142032151e+00, -1.666122363914468e+01 };
        
        double px2,val;
        px2 = 1.0 / ( x * x );
        val = R[5];
        for(int n=4; n>=0; n--) val = val * px2 + R[n];
        val = val / x;
        return val;
    };
    /* */
};

double hyperg_U_CF1(const double a, const int N, const double x ) {
    /* continued fraction for (a+N) U(a+N+1,3/2,x)/U(a+N,3/2,x) */
    /*          some of the code is taken from GSL 1.9          */
    const double RECUR_BIG = 2.0000e+64;
    const int mxt = 1000;
    
    double Anm2,Bnm2,Anm1,Bnm1,An,Bn;
    double a1,b1,fn,old_fn,del,an,bn;
    
    Anm2 = 1.0;
    Bnm2 = 0.0;
    Anm1 = 0.0;
    Bnm1 = 1.0;
    a1   = -( a + N );
    b1   =  ( 1.5 - 2.0 * a - x - 2.0 * ( N + 1 ) );
    An   = b1 * Anm1 + a1 * Anm2;
    Bn   = b1 * Bnm1 + a1 * Bnm2;
    fn   = An / Bn;
    
    /* Steed's algorithm for the continued fraction */
    for(int n=2; n<=mxt; n++) {
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        an  = -( a + N + n - 1.5 ) * ( a + N + n - 1.0 );
        bn  =  ( 1.5 - 2.0 * a - x - 2.0 * ( N + n ) );
        An  = bn * Anm1 + an * Anm2;
        Bn  = bn * Bnm1 + an * Bnm2;
    
        if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
            /* scale An and Bn to avoid an overflow */
            An   /= RECUR_BIG;
            Bn   /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
            Bnm2 /= RECUR_BIG;
        };
    
        old_fn = fn;
        fn     = An / Bn;
        del    = old_fn / fn;
    
        if(fabs(del - 1.0) < numeric_limits<double>::epsilon()) break;
    };
    return fn;
    /* */
};

void SphericalBesselJ( const int nmax, const double z, double* jn ) {
    /* spherical Bessel functions, j_n(z), z>0 */
    if( nmax > 25 ) {
        cout << " The code has never been tested for n>25!" << endl;
        cout << " Emergency halt." << endl;
        exit( EXIT_FAILURE );
    };
    
    if( z < numeric_limits<double>::epsilon() ) {
        jn[0] = 1.0;
        return ;
    };
    
    if( z < 25.0 ) {
        usint n_strt,n;
        n_strt = 51;
        double* tn = new double[n_strt];
        n_strt--;
        
        tn[n_strt] = 0.0;
        for(n=n_strt; n>=1; n--) tn[n-1] = z / ( 2*n + 1 - z * tn[n] );
    
        jn[0] = sin(z) / z;
        for(n=1; n<=nmax; n++) jn[n] = tn[n-1] * jn[n-1];
        
        delete []tn;
    } else {
        usint n;
        jn[0] = sin(z) / z;
        jn[1] = ( jn[0] - cos(z) ) / z;
        for(n=1; n<nmax; n++) jn[n+1] = ( 2*n + 1 ) * jn[n] / z - jn[n-1];
    };
    /* */
};

void BuildAlphaIJ( const usint ij_max, double** v ) {
    /* \alpha_{ij} = \int_0^\pi d\phi \sin^i \phi \cos^j \phi */
    usint i,j;
    
    v[0][0] = M_PI;
    v[1][0] = 2.0;
    for(i=2; i<=ij_max; i++) v[i][0] = v[i-2][0] * ( i - 1 ) / i;
    
    for(j=2; j<=ij_max; j+=2)
        for(i=0; i<=ij_max-j; i++)
            v[i][j] = v[i][j-2] - v[i+2][j-2];
    /* */
};

void CartSolidAng( const usint i, const usint j, const usint k, 
                   const usint l, double** a, double* proj ) {
    /* \int d\Omega x^i y^j z^k Y_{lm}(r) / r^(i+j+k) */
    double sum;
    int m;
    usint x,y,z;
    
    for(m=-l; m<=l; m++) {
        sum = 0.0;
        for(x=0; x<=l; x++)
            for(y=0; y<=l-x; y++) 
                for(z=0; z<=l-x-y; z++)
            sum += NoNormCalcClmR( l, m, x, y, z )
                 * CartSolidAng( i + x, j + y, k + z, a );
        proj[m+l] = sum;
    };
    /* */
};

void CompSolidHarm( const double x, 
                    const double y,
                    const double z,
                    const usint l, double* ylm ) {
    /* real solid harmonics (no Racah) at [x,y,z] */
    double sum;
    double xp,yp,zp,rr;
    double xr,yr,zr;
    int m;
    usint i,j,k;
    
    rr = x * x + y * y + z * z;
    rr = sqrt( rr );
    
    if( rr < numeric_limits<double>::epsilon() ) {
        memset( ylm, 0, (l+l+1) * sizeof(double) );
        if( l == 0 ) ylm[0] = 1.0 / sqrt_4pi; // by convention
        return ;
    };
    
    xr = x / rr;
    yr = y / rr;
    zr = z / rr;
    
    for(m=-l; m<=l; m++) {
        sum = 0.0;
        xp  = 1.0;
        for(i=0; i<=l; i++) {
            yp = 1.0;
            for(j=0; j<=l-i; j++) {
                zp = 1.0;
                for(k=0; k<=l-i-j; k++) {
                    sum += NoNormCalcClmR( l, m, i, j, k )
                         * xp * yp * zp;
                    zp = zp * zr;
                };
                yp = yp * yr;
            };
            xp = xp * xr;
        };
        ylm[m+l] = sum;
    };
    /* */
};

void CompAngFac( const double kx, 
                 const double ky,
                 const double kz,
                 const usint lP, const usint lG, 
                 const usint L , double** aij,
                 double* ang_fac ) {
    /* angular factors for the projector operators */
    int M;
    usint a,b,c;
    usint ia,ja,ka,mi;
    double s1,s2;
    arr1d <double> ylm( 2*L + 1 );
    arr1d <double> ang( 2*L + 1 );
    
    CompSolidHarm( kx, ky, kz, L, ylm.v );
    
    ia = 0;
    ja = 0;
    for(mi=0; mi<crt_siz[lG]; mi++) {
        ka = lG - ia - ja;
        
        s1 = 0.0;
        for(a=0; a<=lP; a++) for(b=0; b<=lP-a; b++) for(c=0; c<=lP-a-b; c++) {
            s2 = 0.0;
            CartSolidAng( ia + a, ja + b, ka + c, L, aij, ang.v );
            for(M=-L; M<=L; M++) s2 += ang[M+L] * ylm[M+L];
            s1 += s2 * NoNormCalcClmR( lP, 0, a, b, c );
        };
        
        ang_fac[mi] = s1;
            
        ja++;
        if( ja > lG - ia ) { ja = 0; ia++; };
    };
    /* */
};

void CAPAngFac( const double kx, const double ky, const double kz,
                const usint lG, const usint Lmax , 
                double** aij  , double** ang_fac ) {
    /* angular factors for the atomic CAP integrals */
    int M;
    usint ia,ja,ka,mi;
    usint L,twoL;
    twoL = 2*Lmax + 1;
    
    arr1d <double> ylm( twoL );
    arr1d <double> ang( twoL );
    
    for(L=0; L<=Lmax; L++) {
        ylm.zero(); ang.zero();
        CompSolidHarm( kx, ky, kz, L, ylm.v );
        
        ia = 0;
        ja = 0;
        for(mi=0; mi<crt_siz[lG]; mi++) {
            ka = lG - ia - ja;
            
            CartSolidAng( ia, ja, ka, L, aij, ang.v );
            for(M=0; M<twoL; M++) 
                ang_fac[L][mi] += ang[M] * ylm[M];
            
            ja++;
            if( ja > lG - ia ) { ja = 0; ia++; };
        };
    };
    /* */
};
