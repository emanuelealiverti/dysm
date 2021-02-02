# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "skew_cdf.hpp"

//****************************************************************************80

double alnorm ( double x, bool upper )

//****************************************************************************80
//
//  Purpose:
//
//    ALNORM computes the cumulative density of the standard normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by David Hill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
{
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  bool up;
  double utzero = 18.66;
  double value;
  double y;
  double z;

  up = upper;
  z = x;

  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y 
      / ( y + a1 + b1 
      / ( y + a2 + b2 
      / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y ) 
      / ( z + c1 + d1 
      / ( z + c2 + d2 
      / ( z + c3 + d3 
      / ( z + c4 + d4 
      / ( z + c5 + d5 
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

  return value;
}



double tfn ( double x, double fx )

//****************************************************************************80
//
//  Purpose:
//
//    TFN calculates the T-function of Owen.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by JC Young, Christoph Minder.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    MA Porter, DJ Winstanley,
//    Remark AS R30:
//    A Remark on Algorithm AS76:
//    An Integral Useful in Calculating Noncentral T and Bivariate
//    Normal Probabilities,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, page 113.
//
//    JC Young, Christoph Minder,
//    Algorithm AS 76: 
//    An Algorithm Useful in Calculating Non-Central T and 
//    Bivariate Normal Distributions,
//    Applied Statistics,
//    Volume 23, Number 3, 1974, pages 455-457.
//
//  Parameters:
//
//    Input, double X, FX, the parameters of the function.
//
//    Output, double TFN, the value of the T-function.
//
{
# define NG 5

  double fxs;
  int i;
  double r[NG] = {
    0.1477621, 
    0.1346334, 
    0.1095432, 
    0.0747257, 
    0.0333357 };
  double r1;
  double r2;
  double rt;
  double tp = 0.159155;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  double tv4 = 1.0E-05;
  double u[NG] = {
    0.0744372, 
    0.2166977, 
    0.3397048,
    0.4325317, 
    0.4869533 };
  double value;
  double x1;
  double x2;
  double xs;
//
//  Test for X near zero.
//
  if ( fabs ( x ) < tv1 )
  {
    value = tp * atan ( fx );
    return value;
  }
//
//  Test for large values of abs(X).
//
  if ( tv2 < fabs ( x ) )
  {
    value = 0.0;
    return value;
  }
//
//  Test for FX near zero.
//
  if ( fabs ( fx ) < tv1 )
  {
    value = 0.0;
    return value;
  }
//
//  Test whether abs ( FX ) is so large that it must be truncated.
//
  xs = - 0.5 * x * x;
  x2 = fx;
  fxs = fx * fx;
//
//  Computation of truncation point by Newton iteration.
//
  if ( tv3 <= log ( 1.0 + fxs ) - xs * fxs )
  {
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;

    for ( ; ; )
    {
      rt = fxs + 1.0;

      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) 
      / ( 2.0 * x1 * ( 1.0 / rt - xs ) );

      fxs = x2 * x2;

      if ( fabs ( x2 - x1 ) < tv4 )
      {
        break;
      }
      x1 = x2;
    }
  }
//
//  Gaussian quadrature.
//
  rt = 0.0;
  for ( i = 0; i < NG; i++ )
  {
    r1 = 1.0 + fxs * pow ( 0.5 + u[i], 2 );
    r2 = 1.0 + fxs * pow ( 0.5 - u[i], 2 );

    rt = rt + r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
  }

  value = rt * x2 * tp;

  return value;
# undef NG
}

void psnc( double *out, double *x, double xi, double omega, double alpha, int n)
{
	double xc;

   for(int i = 0; i < n; i++) {

	xc = (x[i]-xi)/omega;
	out[i] = -2.0*tfn(xc,alpha);
	out[i] += alnorm(xc,false);
   }

	return;
}

void diffc(double *out, double *x, int n) {
	out[0] = x[0];
   for(int i = 1; i < n; i++) {
	   out[i] = x[i] - x[i-1];
   }
   return;
}

