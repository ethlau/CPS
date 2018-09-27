////////////////////////////////////////////////////////////////////////////////
// File: power_series_Si.c                                                    //
// Routine(s):                                                                //
//    xPower_Series_Si                                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xPower_Series_Si( long double x )                              //
//                                                                            //
//  Description:                                                              //
//     The sin integral, Si(x), is the integral with integrand                //
//                             sin(t) / t                                     //
//     where the integral extends from 0 to x.                                //
//     The power series representation for Si(x) is                           //
//               x Sum (-x^2)^j / [(2j+1) (2j+1)!] where the sum extends      //
//     over j = 0, ,,,.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the sin integral Si().                      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the sin integral Si evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xPower_Series_Si( x );                                             //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>           // required for fabsl(), expl()

long double xPower_Series_Si( long double x )
{ 
   long double sum;
   long double xx = - x * x;
   int n = (int) (3.42 * xx + 7.46 * fabsl(x) + 6.95);
   int k = n + n;

       // If the argument is sufficiently small use the approximation //
                         // Si(x) = x exp(-x^2/18) //
  
   if ( fabsl(x) <= 0.0003L ) return x * expl(xx / 18.0L);

        // Otherwise evaluate the power series expansion for Si(x). //

   sum = xx / ( (k + 1) * (k + 1) * k) + 1.0L / (k - 1);
   for (k--; k >= 3; k -=2) {
      sum *= xx / ( k * (k - 1) );
      sum += 1.0L / (k - 2);
   }
   return x * sum;
}
