////////////////////////////////////////////////////////////////////////////////
// File: power_series_Cin.c                                                   //
// Routine(s):                                                                //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xPower_Series_Cin( long double x )                             //
//                                                                            //
//  Description:                                                              //
//     The entire cos integral, Cin(x), is the integral with integrand        //
//                             (1 -  cos(t)) / t                              //
//     where the integral extends from 0 to x.                                //
//     The power series representation for Cin(x) is                          //
//                 - Sum (-x^2)^j / [(2j) (2j)!] where the sum extends        //
//     over j = 1, ,,,.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     long double x   The argument of the entire cos integral Cin().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the entire cos integral Cin evaluated at x.               //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xPower_Series_Cin( x );                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                      // required for fabsl(), expl()

long double xPower_Series_Cin( long double x )
{ 
   long double sum;
   long double xx = - x * x;
   int n = (int) (2.41 * xx + 7.15 * fabsl(x) + 7.00);
   int k = n + n;
  
       // If the argument is sufficiently small use the approximation //
                      // Cin(x) = x^2 exp(-x^2/24) / 4 //

   if ( fabsl(x) < 0.00025L ) return - (expl(xx / 24.0L) * xx) / 4.0L;

       // Otherwise evaluate the power series expansion for Cin(x). //

   sum = xx / ( k * k * (k - 1)) + 1.0L / (k - 2);
   for (k -= 2; k > 2; k -= 2) {
      sum *= xx / ( k * (k - 1) );
      sum += 1.0L / (k - 2);
   }
   return -xx * sum / 2.0L;
}
