
void Faddeeva( const double xi, const double yi,
               double & u, double & v, bool & flag ) {
    /* Faddeeva function for complex arguments;
     * based on: G. M. Poppe and C. M. J. Wijers,
     * ACM Trans. Math. Soft. 16, 38 (1990).
     */
    double factor,rmaxreal,rmaxexp;
    double xabs,yabs,x,y,qrho,xabsq,xquad,yquad;
    double xsum,ysum,xaux,daux,u1,u2,v1,v2;
    double h,h2,qlambda,rx,sx,ry,sy,tx,ty,c,w1;
    bool a,b;
    int n,j,kapn,nu,np1;
    
    factor   = 2.0 / sqrt( M_PI );
    rmaxreal = sqrt( numeric_limits<double>::max() );
    rmaxexp  = 2.0 * log( rmaxreal ) - M_LN2;
    
    flag = false;
    xabs = abs( xi );
    yabs = abs( yi );
    x    = xabs / 6.3;
    y    = yabs / 4.4;
    
    if( xabs > rmaxreal || yabs > rmaxreal ) {
        flag = true; return ;
    };
    
    qrho  = x * x + y * y;
    xabsq = xabs * xabs;
    xquad = xabsq - yabs * yabs;
    yquad = 2.0 * xabs * yabs;
    a     = qrho < 0.085264;
    
    if( a ) {
        qrho = ( 1.0 - 0.85 * y ) * sqrt( qrho );
        n    = (int) nearbyint( 6.0 + 72.0 * qrho );
        j    = 2*n + 1;
        xsum = 1.0 / j;
        ysum = 0.0;
        
        for(int i=n; i>=1; i--) {
            j    = j - 2;
            xaux = ( xsum * xquad - ysum * yquad ) / i;
            ysum = ( xsum * yquad + ysum * xquad ) / i;
            xsum = xaux + 1.0 / j;
        };
        
        u1   = -factor * ( xsum * yabs + ysum * xabs ) + 1.0;
        v1   = +factor * ( xsum * xabs - ysum * yabs );
        daux = +exp( -xquad );
        u2   = +daux * cos( yquad );
        v2   = -daux * sin( yquad );
        
        u    = u1 * u2 - v1 * v2;
        v    = u1 * v2 + v1 * u2;
    } else {
        h2      = 0.0;
        qlambda = 0.0;
        if( qrho > 1.0 ) {
            h    = 0.0;
            kapn = 0;
            qrho = sqrt( qrho );
            nu   = (int) nearbyint( 3.0 + ( 1442.0 / ( 26.0 * qrho + 77.0 ) ) );
        } else {
            qrho = ( 1.0 - y ) * sqrt( 1.0 - qrho );
            h    = 1.88 * qrho;
            h2   = 2.0 * h;
            kapn = (int) nearbyint(  7.0 + 34.0 * qrho );
            nu   = (int) nearbyint( 16.0 + 26.0 * qrho );
        };
        
        b = h > 0.0;
        
        if( b ) qlambda = pow( h2, kapn );
        
        rx = 0.0;
        ry = 0.0;
        sx = 0.0;
        sy = 0.0;
        
        for(n=nu; n>=0; n--) {
            np1 = n + 1;
            tx  = yabs + h + np1 * rx;
            ty  = xabs - np1 * ry;
            c   = 0.5 / ( tx * tx + ty * ty );
            rx  = c * tx;
            ry  = c * ty;
            if( b && ( n <= kapn ) ) {
                tx = qlambda + sx;
                sx = rx * tx - ry * sy;
                sy = ry * tx + rx * sy;
                qlambda = qlambda / h2;
            };
        };
        
        if( abs( h ) < numeric_limits<double>::epsilon() ) {
            u = factor * rx;
            v = factor * ry;
        } else {
            u = factor * sx;
            v = factor * sy;
        };
        
        if( abs( yabs ) < numeric_limits<double>::epsilon() ) 
            u = exp( -xabs * xabs );
    };
    
    if( yi < 0.0 ) {
        if( a ) {
            u2 = 2.0 * u2;
            v2 = 2.0 * v2;
        } else {
            xquad = -xquad;
            if( xquad > rmaxexp ) {
                flag = true; return ;
            };
            w1 = +2.0 * exp( xquad );
            u2 = +w1 * cos( yquad );
            v2 = -w1 * sin( yquad );
        };
        
        u = u2 - u;
        v = v2 - v;
        if( xi > 0.0 ) v = -v;
    } else {
        if( xi < 0.0 ) v = -v;
    };
    /* */
};

void Faddeeva( const quad xi, const quad yi,
               quad & u, quad & v, bool & flag ) {
    /* Faddeeva function for complex arguments,
     * quadruple precision subroutine version;
     * based on: G. M. Poppe and C. M. J. Wijers,
     * ACM Trans. Math. Soft. 16, 38 (1990).
     */
    quad factor,rmaxreal,rmaxexp;
    quad xabs,yabs,x,y,qrho,xabsq,xquad,yquad;
    quad xsum,ysum,xaux,daux,u1,u2,v1,v2;
    quad h,h2,qlambda,rx,sx,ry,sy,tx,ty,c,w1;
    bool a,b;
    int n,j,kapn,nu,np1;
    
    factor   = 2.0q / sqrtq( M_PIq );
    rmaxreal = sqrt( numeric_limits<double>::max() );
    rmaxexp  = 2.0q * logq( rmaxreal ) - M_LN2q;
    
    flag = false;
    xabs = fabsq( xi );
    yabs = fabsq( yi );
    x    = xabs / 6.3q;
    y    = yabs / 4.4q;
    
    if( xabs > rmaxreal || yabs > rmaxreal ) {
        flag = true; return ;
    };
    
    qrho  = x * x + y * y;
    xabsq = xabs * xabs;
    xquad = xabsq - yabs * yabs;
    yquad = 2.0q * xabs * yabs;
    a     = qrho < 0.085264q;
    
    PrintQuadMath( qrho );
    
    if( a ) {
        qrho = ( 1.0q - 0.85q * y ) * sqrtq( qrho );
        n    = (int) nearbyint( 16.0 + 132.0 * (double) qrho );
        j    = 2*n + 1;
        xsum = 1.0q / j;
        ysum = 0.0q;
        
        for(int i=n; i>=1; i--) {
            j    = j - 2;
            xaux = ( xsum * xquad - ysum * yquad ) / i;
            ysum = ( xsum * yquad + ysum * xquad ) / i;
            xsum = xaux + 1.0q / j;
        };
        
        u1   = -factor * ( xsum * yabs + ysum * xabs ) + 1.0q;
        v1   = +factor * ( xsum * xabs - ysum * yabs );
        daux = +expq( -xquad );
        u2   = +daux * cosq( yquad );
        v2   = -daux * sinq( yquad );
        
        u    = u1 * u2 - v1 * v2;
        v    = u1 * v2 + v1 * u2;
    } else {
        h2      = 0.0q;
        qlambda = 0.0q;
        if( qrho > 1.0q ) {
            h    = 0.0q;
            kapn = 0;
            qrho = sqrtq( qrho );
            nu   = (int) nearbyint( 12.0 + ( 2800.0 / ( 26.0 * (double) qrho + 77.0 ) ) );
        } else {
            qrho = ( 1.0q - y ) * sqrtq( 1.0q - qrho );
            h    = 1.88q * qrho;
            h2   = 2.0q * h;
            kapn = (int) nearbyint(  7.0 + 34.0 * (double) qrho );
            nu   = (int) nearbyint( 16.0 + 26.0 * (double) qrho );
        };
        
        b = h > 0.0;
        
        if( b ) qlambda = powq( h2, kapn );
        
        rx = 0.0;
        ry = 0.0;
        sx = 0.0;
        sy = 0.0;
        
        for(n=nu; n>=0; n--) {
            np1 = n + 1;
            tx  = yabs + h + np1 * rx;
            ty  = xabs - np1 * ry;
            c   = 0.5q / ( tx * tx + ty * ty );
            rx  = c * tx;
            ry  = c * ty;
            if( b && ( n <= kapn ) ) {
                tx = qlambda + sx;
                sx = rx * tx - ry * sy;
                sy = ry * tx + rx * sy;
                qlambda = qlambda / h2;
            };
        };
        
        if( fabsq( h ) < numeric_limits<double>::epsilon()
                       * numeric_limits<double>::epsilon() ) {
            u = factor * rx;
            v = factor * ry;
        } else {
            u = factor * sx;
            v = factor * sy;
        };
        
        if( fabsq( yabs ) < numeric_limits<double>::epsilon()
                          * numeric_limits<double>::epsilon() ) 
            u = expq( -xabs * xabs );
    };
    
    if( yi < 0.0q ) {
        if( a ) {
            u2 = 2.0q * u2;
            v2 = 2.0q * v2;
        } else {
            xquad = -xquad;
            if( xquad > rmaxexp ) {
                flag = true; return ;
            };
            w1 = +2.0 * expq( xquad );
            u2 = +w1  * cosq( yquad );
            v2 = -w1  * sinq( yquad );
        };
        
        u = u2 - u;
        v = v2 - v;
        if( xi > 0.0q ) v = -v;
    } else {
        if( xi < 0.0q ) v = -v;
    };
    /* */
};

void CAPRadFac_General( const int lmax, const int nmax,
                const double k, const double p,
                double** rad_fac ) {
    /* radial factors for the atomic CAP integrals */
    if( nmax == 0 ) return ;
    bool flag;
    int l,n,ln;
    double sq_p,j01,j02,j03,fac,last;
    double fed_re,fed_im,arg_re,arg_im;
    cdouble fed,dfed,tmp,face;
    
    /* special case for small k < 1 */
    if( k * k / p < 1.0 ) {
        sq_p   = sqrt( p );
        fed_re = exp( -p );
        fed_im = 2.0 * p;
        arg_re = -k * k / 2.0;
        
        arr1d <double> gn( lmax + nmax + 1 );
        gn[0] = 0.25 * sqrt_4pi * erfc( sq_p ) / sq_p;
        gn[1] = fed_re / fed_im;
        
        for(n=2; n<=nmax+lmax; n++)
            gn[n] = ( ( n - 1 ) * gn[n-2] + fed_re ) / fed_im;
        
        for(l=0; l<=lmax; l++) for(n=l+1; n<=nmax; n++) {
            ln     = l + n;
            last   = 0.0;
            arg_im = 1.0;
            fac    = gn[ln];
            
            for(int m=0; m<=50; m++) {
                j01    = arg_im * fac / dfact[l+m+1] / fact[m];
                last  += j01;
                if( abs( j01 / last ) < numeric_limits<double>::epsilon() ) break;
                fac    = ( ( ln + m + m + 1 ) * fac + fed_re ) / fed_im;
                arg_im = arg_im * arg_re;
            };
            
            rad_fac[l][n] = last * pow( k, l );
        };
        
        return ;
    };
    
    /* part 1: J_{01} and J_{02} from scratch */
    flag   = false;
    sq_p   = sqrt( p );
    arg_re = k / ( 2.0 * sq_p );
    arg_im = sq_p;
    
    Faddeeva( arg_re, arg_im, fed_re, fed_im, flag );
    fed  = fed_re + I * fed_im;
    dfed = 2.0 * ( 2.0 * I / sqrt_4pi - fed * ( arg_re + I * arg_im ) );
    
    face = exp( I * k - p );
    tmp  = fed * face;
    fac  = 0.25 * sqrt_4pi / sq_p / k;
    j01  = fac * tmp.imag();
    
    tmp  = face * ( I * fed + 0.5 / sq_p * dfed );
    j02  = -fac * tmp.real();
    
    /* part 2: recursion to increase n */
    rad_fac[0][0] = 0.0;
    rad_fac[0][1] = j01;
    if( nmax == 1 ) return ;
    
    rad_fac[0][2] = j02;
    if( nmax == 2 ) return ;
    
    fed_re = exp( -p ) / ( 4.0 * p * p );
    fed_im = 1.0 / ( 2.0 * p );
    arg_im = sin( k );
    arg_re = ( cos( k ) + arg_im / k / fed_im ) * fed_re;
    arg_im = -arg_im * fed_re / k;
    fed_re = k * fed_im;
    fed_re = fed_re  * fed_re;
    
    j03 = ( 0.5 / p - fed_re ) * j01 + arg_re;
    rad_fac[0][3] = j03;
    
    for(n=1; n<=nmax-3; n++) 
        rad_fac[0][n+3] = ( fed_im * ( n + n + 1 ) - fed_re ) * rad_fac[0][n+1]
                        - n * ( n - 1 ) * rad_fac[0][n-1] * fed_im * fed_im
                        + arg_re + arg_im * n;
    n = nmax - 2;
    last = ( fed_im * ( n + n + 1 ) - fed_re ) * rad_fac[0][n+1]
         - n * ( n - 1 ) * rad_fac[0][n-1] * fed_im * fed_im
         + arg_re + arg_im * n;
    
    /* part 3: build l at the cost of n */
    fed_re = 2.0 * p;
    fed_im = exp( -p ) * sin( k ) / k;
    
    if( lmax > 0 ) {
        for(n=2; n<nmax; n++)
            rad_fac[1][n] = ( fed_im + n * rad_fac[0][n-1] 
                          -   fed_re * rad_fac[0][n+1] ) / k;
    
        rad_fac[1][nmax]  = ( nmax * rad_fac[0][nmax-1] 
                          -   fed_re * last + fed_im ) / k;
    };
    
    for(l=1; l<lmax; l++) for(n=l+2; n<=nmax; n++)
        rad_fac[l+1][n] = ( l + l + 1 ) * rad_fac[l][n-1] / k - rad_fac[l-1][n];
    
    /* */
};

