#include <math.h>
#include <stdio.h>
#include "proto.h"
#include "allvars.h"
#include <math.h>

double n_ps(double mass, double zred)
{
  //PS mass function
  // dN/dM in a comoving volume (Mpc/h)**3
  
  double sigmam,sigmam_linear,dsdm;
  double delta_c,a_present,a,deltac_z;
  double rho_c=2.7755e11;
  double res,D,a_init;
  
  //convert mass in units of Msun
  
  mass=mass/HubbleParam;
  
  sigmam = sigmam_linear_eh(mass);
  dsdm   = dsdm_eh(mass);
  
  delta_c  = 3.0/20.0*pow(12.0*PI,2.0/3.0);

  a_present=1.0; a_init=1.0e-4; 
  a=1.0/(1.0+zred);
  
  D=GrowthFactor(a_init,a)/GrowthFactor(a_init,a_present);
  
  deltac_z = delta_c/D;
  
  res = sqrt(2.0/PI)*rho_c*Omega*HubbleParam*HubbleParam/mass
    *deltac_z/sigmam/sigmam
    *(-dsdm)*exp(-0.5*deltac_z*deltac_z/sigmam/sigmam);

  
  //from /(Mpc)^3 to /(Mpc/h)*3
  res = res /HubbleParam/HubbleParam/HubbleParam;
  
  //compensate Msun
  res = res /HubbleParam;

  // avoid underflow

  if(res < 1.0e-60)
    res = 1.e-60;
  
  return res;
}

double n_st(double mass, double zred)
{
  //ST mass function
  // dN/dM in a comoving volume (Mpc/h)**3
  
  double sigmam,sigmam_linear,dsdm;
  double delta_c,a_present,a,deltac_z;
  double rho_c=2.7755e11;
  double Ac=0.32222;
  double as=0.707;
  double  p=0.3;
  double res,D,a_init;
  
  //convert mass in units of Msun
  
  mass=mass/HubbleParam;
  
  sigmam = sigmam_linear_eh(mass);
  dsdm   = dsdm_eh(mass);
  
  delta_c  = 3.0/20.0*pow(12.0*PI,2.0/3.0);

  a_present=1.0; a_init=1.0e-5; 
  a=1.0/(1.0+zred);
  
  D=GrowthFactor(a_init,a)/GrowthFactor(a_init,a_present);
  
  deltac_z = delta_c/D;
  
  res = Ac*sqrt(2.0*as/PI)*rho_c*Omega*HubbleParam*HubbleParam/mass
    *(1.0+pow(sigmam*sigmam/deltac_z/deltac_z/as,p))
    *deltac_z/sigmam/sigmam
    *(-dsdm)*exp(-0.5*as*deltac_z*deltac_z/sigmam/sigmam);

  
  //from /(Mpc)^3 to /(Mpc/h)*3
  res = res /HubbleParam/HubbleParam/HubbleParam;
  
  //compensate Msun
  res = res /HubbleParam;

  // avoid underflow

  if(res < 1.0e-60)
    res = 1.e-60;
  
  return res;
}


double n_yahagi(double mass, double zred)
{
  // Yahagi et al. (2004) mass function; eq (7)
  // dN/dM in a comoving volume (Mpc/h)**3
  
  double sigmam,sigmam_linear,dsdm;
  double delta_c,a_present,a,deltac_z;
  double rho_c=2.7755e11;
  double Ac=0.298;
  double Bc=0.893;
  double Cc=1.39;
  double Dc=0.408;
  double res,D,a_init;
  
  //convert mass in units of Msun
  
  mass=mass/HubbleParam;
  
  sigmam = sigmam_linear_eh(mass);
  dsdm   = dsdm_eh(mass);
  
  delta_c  = 3.0/20.0*pow(12.0*PI,2.0/3.0);

  a_present=1.0; a_init=1.0e-5; 
  a=1.0/(1.0+zred);
  
  D=GrowthFactor(a_init,a)/GrowthFactor(a_init,a_present);
  
  deltac_z = delta_c/D;
  
  res = Ac*rho_c*Omega*HubbleParam*HubbleParam/mass
    *(1.0+pow((Bc*deltac_z/sigmam/sqrt(2.0)),Cc))
    *pow((deltac_z/sigmam), Dc)/sigmam
    *(-dsdm)*exp(-0.5*Bc*Bc*deltac_z*deltac_z/sigmam/sigmam);

  
  //from /(Mpc)^3 to /(Mpc/h)*3
  res = res /HubbleParam/HubbleParam/HubbleParam;
  
  //compensate Msun
  res = res /HubbleParam;

  // avoid underflow

  if(res < 1.0e-60)
    res = 1.e-60;
  
  return res;
}



double n_jk(double mass, double zred)
{
  //PS mass function
  // dN/dM in a comoving volume (Mpc/h)**3
  
  double sigmam,sigmam_linear,dsdm;
  double delta_c,a_present,a,deltac_z;
  double rho_c=2.7755e11;
  double Ac=0.32222;
  double as=0.707;
  double  p=0.3;
  double res,D,a_init;
  
  //convert mass in units of Msun
  
  mass=mass/HubbleParam;
  
  sigmam = sigmam_linear_eh(mass);
  dsdm   = dsdm_eh(mass);
  
  delta_c  = 3.0/20.0*pow(12.0*PI,2.0/3.0);

  a_present=1.0; a_init=1.0e-5; 
  a=1.0/(1.0+zred);
  
  D=GrowthFactor(a_init,a)/GrowthFactor(a_init,a_present);
  
  deltac_z = delta_c/D;
  
  //  res = 0.315*rho_c*Omega*HubbleParam*HubbleParam/mass
  //    *1.0/sigmam
  //    *(-dsdm)*exp(-pow(fabs(log(1.0/sigmam/D)+0.61),3.8));

  //Jenkins (2001) equation B2
  res = 0.301*rho_c*Omega*HubbleParam*HubbleParam/mass
    *1.0/sigmam
    *(-dsdm)*exp(-pow(fabs(log(1.0/sigmam/D)+0.64),3.88));

  
  //from /(Mpc)^3 to /(Mpc/h)*3
  res = res /HubbleParam/HubbleParam/HubbleParam;
  
  //compensate Msun
  res = res /HubbleParam;

  // avoid underflow

  if(res < 1.0e-60)
    res = 1.e-60;
  
  return res;
}






double sigmam_ks(double mass)
{ //      sigma(M) at z=0 extrapolated by linear theory 

  double mass12,mass_8,sigmam_8,res;
  double  rho_c=2.7755e11, p=0.0873;
    
  //
  // ==================================================================
  //
  // Kitayama and Suto (ApJ 1996, 469, 480, eq.B5)
  // mass   : in units of M_sun
  // mass12 : in units of 10**12 M_sun
  // ----------------------------------------------
  mass=mass/HubbleParam; //Kitayama Suto use Msun

  mass12 = mass*Gamma*Gamma*HubbleParam*HubbleParam/1.e12;
  res = 1.0 + 2.208*pow(mass12,p) 
    - 0.7668*pow(mass12,2*p)
    + 0.7949*pow(mass12,3*p);
    
  res = pow(res,-2.0/9.0/p);
  //cout<<"Sigma8 "<<Sigma8<<" res "<<res<<" mass12 "<<mass12<<endl;

  // mass_8: mass within 8 Mpc/h sphere in units of 10**12 M_sun
  // -------------------------------------------------------------
  mass_8   = 4.0*PI/3.0*pow(8.0/HubbleParam,3)
        *rho_c*Omega*HubbleParam*HubbleParam*Gamma*Gamma*HubbleParam*HubbleParam/1.e12;
  sigmam_8 = 1.0 + 2.208*pow(mass_8,p) 
    - 0.7668*pow(mass_8,2*p)+ 0.7949*pow(mass_8,3*p);
  
  sigmam_8 = pow(sigmam_8,-2.0/9.0/p);

  res = res/sigmam_8*Sigma8;

  return res;
    
}




double dsdm_eh(double mass)   // mass   : in units of M_sun
{
  //Modified Kitayama and Suto (ApJ 1996, 469, 480, eq.B5)
  //mass12 : in units of 10**12 M_sun
  //mass_8: mass within 8 Mpc/h sphere in units of 10**12 M_sun
  
  double mass12, mass_8,res,sigmam_8;
  double rho_c=2.7755e11;
  double p=0.0843;
 
  mass12   = mass*Gamma*Gamma*HubbleParam*HubbleParam/1.e12;
  mass_8   = 4.0*PI/3.0*pow(8.0/HubbleParam,3)
    *rho_c*Omega*HubbleParam*HubbleParam*pow(Gamma*HubbleParam,2)/1.e12;
  
  sigmam_8 = 1.0 + 2.262*pow(mass_8,p)- 0.7828*pow(mass_8,2*p)
           + 0.7349*pow(mass_8,3*p);

  sigmam_8 = pow(sigmam_8,-2.0/9.0/p);

  res = -2.0/9.0/p  * pow((1.0 + 2.262*pow(mass12,p) - 0.7828*pow(mass12,2*p) 
			       + 0.7349*pow(mass12,3.*p)),-2.0/9.0/p-1.0)
    *(2.262*p*pow(mass12,p-1.0)-0.7828*2.0*p*pow(mass12,2.0*p-1.0)
      + 0.7349*3.0*p*pow(mass12,3.0*p-1.0))*pow(Gamma*HubbleParam,2)/1.e12;
    
  res = res/sigmam_8*Sigma8;

  return res;
}

double sigmam_linear_eh(double mass)   // mass   : in units of M_sun
{
  double mass12,mass_8,sigmam_8,res;
  double rho_c=2.7755e11;
  double p=0.0843;

  // mass12 : in units of 10**12 M_sun

  //convert mass in units of M_sum
  
  mass12 = mass*Gamma*Gamma*HubbleParam*HubbleParam/1.e12;

  res = 1.0 + 2.262*pow(mass12,p) - 0.7828*pow(mass12,2.*p)
    + 0.7349*pow(mass12,3.*p);
    
  res = pow(res,-2.0/9.0/p);
  //
  // mass_8: mass within 8 Mpc/h sphere in units of 10**12 M_sun
  
  mass_8   = 4.0*PI/3.0*pow(8.0/HubbleParam,3)
  *rho_c*Omega*HubbleParam*HubbleParam*Gamma*Gamma*HubbleParam*HubbleParam/1.e12;
  
  sigmam_8 = 1.0 + 2.262*pow(mass_8,p) - 0.7828*pow(mass_8,2.*p)
    + 0.7349*pow(mass_8,3.*p);
  
  sigmam_8 = pow(sigmam_8,-2.0/9.0/p);
    
  res = res/sigmam_8*Sigma8;
    
  return res;
}
