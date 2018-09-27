////////////////////////////////////////////////////////////////////////////////
// File: sin_cos_integrals.c                                                  //
// Routine(s):                                                                //
//    Sin_Integral_Si                                                         //
//    xSin_Integral_Si                                                        //
//    Entire_Cos_Integral_Cin                                                 //
//    xEntire_Cos_Integral_Cin                                                //
//    Cos_Integral_Ci                                                         //
//    xCos_Integral_Ci                                                        //
//    Sin_Cos_Integrals_Si_Ci                                                 //
//    xSin_Cos_Integrals_Si_Ci                                                //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>        // required for fabsl(), sinl(), cosl(), and logl()
#include <float.h>       // required for LDBL_EPSILON
#include "proto.h"
#include "allvars.h"

////////////////////////////////////////////////////////////////////////////////
// double Sin_Integral_Si( double x )                                         //
//                                                                            //
//  Description:                                                              //
//     The sin integral, Si(x), is the integral with integrand                //
//                             sin(t) / t                                     //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the sin integral Si().                      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the sin integral Si evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Sin_Integral_Si( x );                                              //
////////////////////////////////////////////////////////////////////////////////
double Sin_Integral_Si( double x )
{
   return (double) xSin_Integral_Si( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xSin_Integral_Si( long double x )                              //
//                                                                            //
//  Description:                                                              //
//     The sin integral, Si(x), is the integral with integrand                //
//                             sin(t) / t                                     //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the sin integral Si().                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the sin integral Si evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xSin_Integral_Si( x );                                             //
////////////////////////////////////////////////////////////////////////////////

long double xSin_Integral_Si( long double x )
{
   long double fi, gi, si;

   if ( fabsl(x) <= 1.0L ) return xPower_Series_Si(x);
   xAuxiliary_Sin_Cos_Integrals_fi_gi(fabsl(x), &fi, &gi);
   si = pi2 - cosl(x) * fi - sinl(fabs(x)) * gi;
   return (x < 0.0L) ? -si : si;
}


////////////////////////////////////////////////////////////////////////////////
// double Entire_Cos_Integral_Cin( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The entire cos integral, Cin(x), is the integral with integrand        //
//                             (1 -  cos(t)) / t                              //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the entire cos integral Cin().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the entire cos integral Cin() evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Entire_Cos_Integral_Cin( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Entire_Cos_Integral_Cin( double x )
{
   return (double) xEntire_Cos_Integral_Cin( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xEntire_Cos_Integral_Cin( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The entire cos integral, Cin(x), is the integral with integrand        //
//                             (1 -  cos(t)) / t                              //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the entire cos integral Cin().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the entire cos integral Cin() evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xEntire_Cos_Integral_Cin( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xEntire_Cos_Integral_Cin( long double x )
{
   if ( fabsl(x) <= 1.0L ) return xPower_Series_Cin(x);
   return logl(fabsl(x)) + euler_gamma - Asymptotic_Series_Ci(x);
}


////////////////////////////////////////////////////////////////////////////////
// double Cos_Integral_Ci( double x )                                         //
//                                                                            //
//  Description:                                                              //
//     The cos integral, Ci(x), is the integral with integrand                //
//                               - cos(t) / t                                 //
//     where the integral extends from x to inf.                              //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the cos integral Ci().                      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cos integral Ci() evaluated at x.                     //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Cos_Integral_Ci( x );                                              //
////////////////////////////////////////////////////////////////////////////////
double Cos_Integral_Ci( double x )
{
   long double ci = xCos_Integral_Ci( (long double) x);
      
   return (fabsl(ci) < DBL_MAX) ? (double) ci :
                                             (ci < 0.0L) ? -DBL_MAX : DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xCos_Integral_Ci( long double x )                              //
//                                                                            //
//  Description:                                                              //
//     The cos integral, Ci(x), is the integral with integrand                //
//                               - cos(t) / t                                 //
//     where the integral extends from x to inf.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the cos integral Ci().                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cos integral Ci evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xCos_Integral_Ci( x );                                             //
////////////////////////////////////////////////////////////////////////////////

long double xCos_Integral_Ci( long double x )
{
   if (x == 0.0L) return -LDBL_MAX;
   if (fabsl(x) <= 1.0L) 
      return logl(fabsl(x)) + euler_gamma - xPower_Series_Cin(x);
   return Asymptotic_Series_Ci(x);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Asymptotic_Series_Ci( long double x )                   //
//                                                                            //
//  Description:                                                              //
//     For a large argument x, the cosine integral, Ci(x), can be expressed as//
//                    Ci(x) = sin(x) fi(x) - cos(x) gi(x)                     //
//     where fi(x) and gi(x) are the auxiliary sine and cosine integrals which//
//     in turn are approximated using their respective asymptotic expansions. //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the cosine integral Ci().                   //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cosine integral Ci evaluated at x.                    //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Asymptotic_Series_Ci( x );                                         //
////////////////////////////////////////////////////////////////////////////////

static long double Asymptotic_Series_Ci( long double x )
{
   long double fi, gi;

   xAuxiliary_Sin_Cos_Integrals_fi_gi(fabsl(x), &fi, &gi);

   return sinl(fabsl(x)) * fi - cosl(x) * gi;
}


////////////////////////////////////////////////////////////////////////////////
// void Sin_Cos_Integrals_Si_Ci( double x, double *Si, double *Ci )           //
//                                                                            //
//  Description:                                                              //
//     This routine returns both the sin integral Si(x) and the cos integral  //
//     Ci(x) via the addresses in the argument list.  For x = 0, *Ci is set   //
//     - DBL_MAX.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//       The argument of the sin integral Si() and the cos integral Ci().     //
//     double *Si                                                             //
//        The address of the sin integral Si() evaluated at x.                //
//     double *Ci                                                             //
//        The address of the cos integral Ci() evaluated at x.                //
//                                                                            //
//  Return Value:                                                             //
//     Type void.  The results are returned via the addresses in the argument //
//     list.                                                                  //
//                                                                            //
//  Example:                                                                  //
//     double x, Si, Ci;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     Sin_Cos_Integrals_Si_Ci( x, &Si, &Ci );                                //
////////////////////////////////////////////////////////////////////////////////
void Sin_Cos_Integrals_Si_Ci( double x, double *Si, double *Ci )
{
   long double xSi, xCi;

   xSin_Cos_Integrals_Si_Ci( (long double) x, &xSi, &xCi);
   *Si = (double) xSi;
   *Ci = (fabsl(xCi) < DBL_MAX) ? (double) xCi : (xCi < 0) ? -DBL_MAX:DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// void xSin_Cos_Integrals_Si_Ci( long double x, long double *Si,             //
//                                                          long double *Ci ) //
//                                                                            //
//  Description:                                                              //
//     This routine returns both the sin integral Si(x) and the cos integral  //
//     Ci(x) via the addresses in the argument list.  For x = 0, *Ci is set   //
//     - LDBL_MAX.                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//       The argument of the sin integral Si() and the cos integral Ci().     //
//     long double *Si                                                        //
//        The address of the sin integral Si() evaluated at x.                //
//     long double *Ci                                                        //
//        The address of the cos integral Ci() evaluated at x.                //
//                                                                            //
//  Return Value:                                                             //
//     Type void.  The results are returned via the addresses in the argument //
//     list.                                                                  //
//                                                                            //
//  Example:                                                                  //
//     long double x, Si, Ci;                                                 //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     xSin_Cos_Integrals_Si_Ci( x, &Si, &Ci );                               //
////////////////////////////////////////////////////////////////////////////////
void xSin_Cos_Integrals_Si_Ci(long double x, long double *Si, long double *Ci)
{
   long double fi, gi, sx, cx;

   if ( x == 0.0L) {
      *Si = 0.0L;
      *Ci = -LDBL_MAX;
      return;
   }
   if ( fabsl(x) <= 1.0L ) {
      *Si = xPower_Series_Si(x);
      *Ci = logl(fabsl(x)) + euler_gamma - xPower_Series_Cin(x);
      return;
   }
   xAuxiliary_Sin_Cos_Integrals_fi_gi(fabsl(x), &fi, &gi);
   sx = sinl(fabsl(x));
   cx = cosl(x);
   *Si = pi2 - cx * fi - sx * gi;
   *Ci = sx * fi - cx * gi;
   if (x < 0.0L) *Si = - *Si;
}
