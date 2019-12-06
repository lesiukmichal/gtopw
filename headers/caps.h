#ifndef CAPS_H
#define CAPS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

#include "../headers/gtopw.h"

template <class type> 
    class FnInts {
    private:
    static const int alloc = 100;
    public:
    int nmax;
    type a;
    type b;
    type* fn;
    
    FnInts( const int nmax_, const type a_, const type b_ ) :
        nmax(nmax_), a(a_), b(b_) {
        fn = new type[alloc]; memset(fn,0,alloc*sizeof(type));
    };
    
   ~FnInts() { delete []fn; };
    void clear() { memset(fn,0,alloc*sizeof(type)); };
};

void FnIntsCalc( FnInts <double> & v );

double hyperg_U_CF1(const double a, const int N, const double x );
double eerfc( const double x );

void SphericalBesselJ( const int nmax, const double z, double* jn );
void BuildAlphaIJ( const usint ij_max, double** v );

inline double 
    CartSolidAng( const usint i, const usint j, 
                  const usint k, double** a ) {
    /* x^i y^j z^k integrated over the solid angle */
    if( ( i + j ) % 2 == 1 ) return 0.0;
    return 2.0 * a[i+j+1][k] * a[j][i];
};

void CartSolidAng( const usint i, const usint j, const usint k, 
                   const usint l, double** a, double* proj );

void CompSolidHarm( const double x, 
                    const double y,
                    const double z,
                    const usint l, double* ylm );

void CompAngFac( const double kx, 
                 const double ky,
                 const double kz,
                 const usint lP, const usint lG, 
                 const usint L , double** aij,
                 double* ang_fac );

void CAPAngFac( const double kx, const double ky, const double kz,
                const usint lG, const usint Lmax , 
                double** aij  , double** ang_fac );

void CAPRadFac( const int lmax, const int nmax,
                const int mpar,
                const double k, const double p,
                double** rad_fac );

static const double u0[51] = {
    +3.9215686274509803921568627450980e-02,-7.8040220447416064262062416241581e-03,
    +7.8035345363307908391782952587178e-04,-3.9405212317041278325085086292627e-05, 
    +3.9969141702990316976339472728698e-08,+9.9363752367857347859529336538759e-08,
    -4.0298990356837854707192811470357e-10,-4.9998789428487696734727264757020e-10, 
    +4.0575941890066591342373252186414e-12,+3.1346129689920464534645546865212e-12,
    -4.0783488516140417319755491654108e-14,-2.1915793292560383804057274373263e-14, 
    +4.0905157741083082211278866462749e-16,+1.6328868829636896307916377356814e-16,
    -4.0924778818340459676153497039982e-18,-1.2662634268456445257414125764597e-18, 
    +4.0826791034289915245396369083227e-20,+1.0075684359591647461024647965406e-20,
    -4.0596587480158560940357619778962e-22,-8.1476752823859289013513356308051e-23, 
    +4.0220885705337538991385897701324e-24,+6.6480358792926857313395969102844e-25,
    -3.9688114402237957412434858450007e-26,-5.4412623599779296519154071134286e-27, 
    +3.8988803183557742252808978224990e-28,+4.4434632707189326402286920480646e-29,
    -3.8115960731676191016749731801367e-30,-3.6005976771587851286013740319728e-31, 
    +3.7065425283541329474274763739952e-32,+2.8767878743964993521908234771242e-33,
    -3.5836170712884552431410141423840e-34,-2.2476692273865262962922364138663e-35, 
    +3.4430551516146793861048279420177e-36,+1.6964150856281372219419858411815e-37,
    -3.2854470899079303036402259435107e-38,-1.2112558588128027296298138507201e-39, 
    +3.1117457954312537225676591802724e-40,+7.8386957431660690584330675133403e-42,
    -2.9232642619675335692118040905959e-42,-4.0830063314867779350323218719595e-44, 
    +2.7216620655712112142028429354235e-44,+8.0210319608030310752081564278726e-47,
    -2.5089205158328067991292049301374e-46,+2.0365690535357915160079890463773e-48, 
    +2.2873065946942749940050545079847e-48,-4.4586158771528262200719351547519e-50,
    -2.0593263303801314890464005831596e-50,+6.4855545496740692316715645207021e-52, 
    +1.8276687708153011153701653995660e-52,-8.1369307686420649986048366726631e-54,
    -1.5951422107019195262008676111737e-54 };

#endif //CAPS_H
