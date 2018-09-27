#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "proto.h"
#include "allvars.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

static char TransTable[100];
static double xp[NPOINTS],yp[NPOINTS],yp2[NPOINTS];

static int WhichSpectrum,np,WhichWindow, OPT,WhichGF;
static int OPT_correc;

static double r_tophat,Rsmooth,Delta_c,fdelta_c,Const_MF;

static double AA,BB,CC;
static double B1,B2,B3,B4,B5;
static double nu,sigma, Omega_z;

static double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
static double G, Hubble;

static double Norm, InitTime, Redshift, a_exp;
static double Dplus; /* growth factor */

static double calc_z;

double Omega_de(double a);
double coeff1(double a);
double coeff2(double a);//u"+coeff1(a)u'+coeff2(a)u=0, u=D/a, '=d/dlna
double RungeKutta(double a_in,double a); //calculate linear density growth eq.

double correction_fac(double logm);
double alpha(double k);
double skewness(double r); //cf. arxiv:1009.5085v1
double skewint(double *k,size_t dim,void *params);
double dskewdlogr(double r);
double skewint2(double *k,size_t dim,void *params);
double dWindow(double x);

int main(int argc, char **argv)
{
  int i,j,ip,iM,ir, idummy, nsnap;
  double dk,k,kmin,kmax,x;
  double sigma_minus1,sigma_0,sigma_1,sigma_p,sigma_v;
  double sigma2_minus1,sigma2_0,sigma2_1,sigma2_p,sigma2_v;
  double MeanDensity,M,R_m,R_mGauss;
  double Mplus,sigma_M_z,R_mplus,lnsigma_inv_plus,dlnsigma_inv;
  double dlogM_over_dlnsigma_inv,a_end;
  double answer,lnsigma_inv,singma_M_zplus,res,ans,fract;
  double TableMass[250],TabledN[250];
  double b_ling;
  char   outname[250];

  FILE *fred;
  FILE *fd;
  FILE *fdt;


  //input parameters======================================
  WhichSpectrum=3;  // 1 Eisenstein-Hu, 2 Bond-Efstathiou, 3 CMBFAST table, 4 BBKS
  WhichGF = 1; // 1 apploximation formula, 2 caluculate linear density growth eq.	
  OPT = 1;          // 0 JK  1 ST  2 PS 3 Yahagi 4 Reed et al (astro-ph/0607150v4) 5 Bhattacharya et al (astro-ph/1005.2239)
  OPT_correc =1; // options for fnl correction 1 LV 2 MVJ

  if(argc!=3){
    fprintf(stderr,"usage:\n > massfunc outputname redshift\n");
    exit(1);
  }
  
  calc_z=atof(argv[2]);

  InitTime=1.0/(1.0+calc_z); // expansion parameter at which MF is computed

  sprintf(outname,argv[1]);
	
  Redshift=1.0/InitTime -1.0;


  cout<<InitTime<<" "<<Redshift<<endl;

  AA=6.4/Gamma*1000.0; 
  BB=3.0/Gamma*1000.0; 
  CC=1.7/Gamma*1000.0;  
  nu=1.13;

  i=set_units();
  i=initialize_powerspectrum(WhichSpectrum);

  
//  printf("Dumping the given initial power spectrum with a factor %g \n",Norm);
//  fdt=fopen("pow.dat","w");
//  for(i=1; i<500; i++){
//    k=pow(10.0,-6.0+0.02*i);
//   
//    x=Dplus * 4*PI*k*k*k*PowerSpec(k);  /* 4 pi k3 P(k) */
//    
//    fprintf(fdt,"%g %g \n",1000.0*k,x/1000.0);
//  }
//  fclose(fdt);


  

  Delta_c=1.68647;
  a_exp=1.0/(1.0+Redshift);

  Omega_z = Omega*pow(1.0+Redshift,3.0)/(Omega*pow(1.0+Redshift,3.0) + OmegaLambda);
  
  fdelta_c=delta_c_func();

  double logm, a1, a2, dlogm = 0.05;

  b_ling=GrowthFactor(InitTime,a_exp)/GrowthFactor(InitTime,1.0);
  Rsmooth=8.0; //here use Mpc/h
  Const_MF = Sigma8 * Sigma8/unnsigma(Rsmooth)*b_ling*b_ling;


  cout<<"Linear Growth Factor="<<b_ling<<endl;
  cout<<"delta_c = "<<1.68647*fdelta_c<<endl;
  cout<<"Const_MF= "<<Const_MF<<" "<<unnsigma(Rsmooth)<<endl;

  fd=fopen(outname,"w");
  for(iM=0;iM< 120;iM++){
    logm = 10.0 + iM * dlogm;
    //logm = 14.0 + iM * dlogm;

	res = dndlogm(logm)*correction_fac(logm);
		
    TableMass[iM]=pow(10.0, logm);
	TabledN[iM]=res;
	printf("%g %g %g\n",TableMass[iM],TabledN[iM],correction_fac(logm));
	fprintf(fd,"%g %g %g\n",TableMass[iM],TabledN[iM],correction_fac(logm));
  }

  fclose(fd);



}





//==================================================
int initialize_powerspectrum(int Spectrum)
{
  double res;
  int i;

  for(i=0;i<NPOINTS;i++)
    xp[i]=yp[i]=yp2[i]=0;

  if(WhichSpectrum==3){
    fprintf(stdout,"initialising...\n");
    readCMB_and_do_spline();
  }
  
  InitTime=1/(1+Redshift);

  AA=6.4/Gamma*1000.0; 
  BB=3.0/Gamma*1000.0; 
  CC=1.7/Gamma*1000.0;  
  nu=1.13;


  B1=2.34/Gamma*1000.0; 
  B2=3.89/Gamma*1000.0; 
  B3=16.1/Gamma*1000.0; 
  B4=5.46/Gamma*1000.0; 
  B5=6.71/Gamma*1000.0; 

  
  Norm=1.0;
  res= TopHatSigma2(8.0);
  Norm=Sigma8*Sigma8/res;
  fprintf(stdout,"Computed Norm=%g \n",Norm);
  //fprintf(stdout,"A_s =%g \n",Norm*(1.5*Omega*10000/C/C)*(1.5*Omega*10000/C/C)/(2*PI*PI));

  Dplus= GrowthFactor(InitTime, 1.0);

  return i;
}




//==================================================
int set_units()    /* ... set some units */
{
  UnitLength_in_cm= 3.085678e24; /* 1.0 Mpc */
  UnitMass_in_g=    1.989e43;    /* 1.0e10 solar masses */ 
  UnitVelocity_in_cm_per_s=1e5;  /* 1 km/sec */

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;

  G=GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);
  Hubble = HUBBLE * UnitTime_in_s;
  
  //  cout<<"UnitLength_in_cm        ="<<UnitLength_in_cm<<endl;
  //  cout<<"UnitMass_in_g           ="<<UnitMass_in_g<<endl;
  //  cout<<"UnitVelocity_in_cm_per_s="<<UnitVelocity_in_cm_per_s<<endl;
  //  cout<<"UnitTime_in_s           ="<<UnitTime_in_s<<endl;
  //  cout<<" "<<endl;

  return 1;
}



//==================================================
double PowerSpec(double k)
{
  switch(WhichSpectrum)
    {
    case 1:
      return PowerSpec_EH(k);
      break;
    case 2:
      return PowerSpec_Efstathiou(k);
      break;
    case 3:
      return PowerSpec_CMBFAST(k);
      break;
    case 4:
      return PowerSpec_BBKS(k);
      break;
    default:
      fprintf(stdout,"Not supported\n");  
    }
}





//==================================================
inline double PowerSpec_Efstathiou(double k)
{
  return Norm*k / pow(1+pow(AA*k+pow(BB*k,1.5)+CC*CC*k*k,nu),2/nu);
}




//==================================================
inline double PowerSpec_BBKS(double k)
{
  return Norm*k * pow(log(1.0+B1*k)/(B1*k),2)/ 
    pow(1+ B2*k + B3*B3*k*k + pow(B4*k,3) + pow(B5*k,4),0.5);
}





//==================================================
double   tk_eh(double k)  
{
  double   q,theta,ommh2,a,s,gamma,L0,C0;
  double   tmp;
  double   omegam, ombh2, hubble;

  /* other input parameters */
  hubble= HubbleParam;

  omegam= Omega;
  ombh2=  OmegaBaryon*HubbleParam*HubbleParam;

  //k*= 1000.0;    /* convert to h/Mpc */
  /*k*= HubbleParam;*/  /* convert to 1/Mpc */


  theta = 2.728/2.7;
  ommh2 = omegam*hubble*hubble;
  s     = 44.5*log(9.83/ommh2)/sqrt( 1.+10.*exp(0.75*log(ombh2)) )*hubble;      
  a     = 1.-0.328*log(431.*ommh2)*ombh2/ommh2
            +0.380*log(22.3*ommh2)*(ombh2/ommh2)*(ombh2/ommh2);
  gamma = a+(1.-a)/(1.+exp(4*log(0.43*k*s)) );
  gamma*= omegam*hubble;
  q     = k*theta*theta/gamma;
  L0    = log( 2.*exp(1.)+1.8*q );
  C0    = 14.2 + 731./(1.+62.5*q);
  tmp   = L0/(L0+C0*q*q);
  return(tmp);
}

//==================================================
inline double PowerSpec_EH(double k)
{
  return Norm*k*pow( tk_eh(k), 2);
}



//==================================================
inline double PowerSpec_CMBFAST(double k)
{
  return Norm*k*pow(transfunc_cmbfast(k), 2);
}


//==================================================
double transfunc_cmbfast(double k)
{
  int i;
  double lk;
  double pow_index;
  
  //k *= 1000.0; /* convert to h/Mpc */
  
  lk=log10(k);
  
  if(lk < xp[0])
    return yp[0]; /* usually should not be executed */

  splint(xp-1, yp-1, yp2-1, np, lk, &pow_index);
  
  return pow(10.0,pow_index);
}

//==================================================
void readCMB_and_do_spline()
{
  int i,iline;
  double yp1,ypn;
  char tfunc_table[80];
  int  errorFlag=0;
  FILE *fd;

  sprintf(tfunc_table,"transfer_nishimichi.dat"); 
  /* sprintf(tfunc_table,"TF.lcdm.total"); */
  
  /* sprintf(tfunc_table,"%s",TransTable); */
  fprintf(stdout,"Reading %s .\n",tfunc_table);
  iline=0;
  if(fd=fopen(tfunc_table,"r")){
    while(!feof(fd)){
      fscanf(fd,"%lf %lf \n",&xp[iline],&yp[iline]);
      fprintf(stdout,"%g %g \n",xp[iline],yp[iline]);
      iline++; 
    }
    fclose(fd);
  }
  else{
    fprintf(stdout,"transfer function file %s not found.\n",tfunc_table);
    errorFlag=1;
  }

  fprintf(stdout,"read in %d data points \n",iline);
  
  np=iline;

  for(i=0;i<iline;i++){
    xp[i]=log10(xp[i]);
    yp[i]=log10(yp[i]);
  }
  
  yp1 = 1.e31;
  ypn = 1.e31;
  
  spline(xp-1, yp-1, iline, yp1, ypn, yp2-1);

  for(i=0;i<iline;i++)
    fprintf(stdout,"%g %g \n",xp[i],yp[i]);

}

//==================================================
inline double TopHatSigma2(double R)
{
  r_tophat= R;

  return qromb(sigma2_int, 0, 1000/R);
}




//==================================================
double sigma2_int(double k)
{
  double kr,kr3,kr2,w,x;
  
  kr=r_tophat*k;
  kr2=kr*kr;
  kr3=kr2*kr;

  if(kr<1e-8) return 0;

  w=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
  x=4*PI*k*k*w*w*PowerSpec(k);
 
  return x;
}




//==================================================
double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}


inline double Hubble_a(double a)
{
  double res;

	res= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda/pow(a,3.0*(1.0+w)));

  return res;
}

double Omega_de(double a){
	double res;
	res = OmegaLambda/(Omega*pow(a,3*w)+OmegaLambda);
	return res;
}	

double coeff1(double a){
	double res;
	res = 0.5*(5-3*w*Omega_de(a));
	return res;
}

double coeff2(double a){
	double res;
	res = 1.5*(1-w)*Omega_de(a);
	return res;
}

double growth(double a)
{
  double hubble_a;

	if(w != -1){
		hubble_a= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda/pow(a,3.0*(1.0+w)));
	}else{
		hubble_a= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda);
	}
	
	switch(WhichGF){
	  case 1:
		return hubble_a*qromb(growth_int, 0, a);
		break;
	  case 2:
		return RungeKutta(1.0/(1.0+50.0),a);
		break;
	  default:
		fprintf(stdout,"Not supported\n"); 
	}
}


inline double growth_int(double a)
{
  return pow(a / (Omega + (1-Omega-OmegaLambda)*a + OmegaLambda/pow(a,3*w)), 1.5);
}



//==================================================
inline double integrand_minus1(double k)
{
  double res, factor;
  
  factor=4.0*PI;
  res=factor*PowerSpec(k)*Window(k*Rsmooth)*Window(k*Rsmooth);
  
  return res;
}


inline double integrand_0(double k)
{
  double res, factor;
  
  factor=4.0*PI;
  res=factor*PowerSpec(k)*Window(k*Rsmooth)*Window(k*Rsmooth)*pow(k,2);
  
  return res;
}

inline double integrand_1(double k)
{
  double res, factor;
  
  factor=4.0*PI;
  res=factor*PowerSpec(k)*Window(k*Rsmooth)*Window(k*Rsmooth)*pow(k,4);
  
  return res;
}

inline double integrand_2(double k)
{
  double res, factor;
  
  factor=4.0*PI;
  res=factor*PowerSpec(k)*Window(k*Rsmooth)*Window(k*Rsmooth)*pow(k,6);
  
  return res;
}

inline double integrand_3(double k)
{
  double res, factor;
  
  factor=4.0*PI;
  res=factor*PowerSpec(k)*Window(k*Rsmooth)*Window(k*Rsmooth)*pow(k,8);
  
  return res;
}


inline double Window(double x)
{
  double res;

  if(WhichWindow==1)
    res=WTopHat(x);
  else
    res=WGaussian(x);

  return res;
}

  

inline double WTopHat(double x)
{
  double window;
	if(sin(x)/x-cos(x) < TINY){
		window =0;
	}else{
		window=3.0*(sin(x)-x*cos(x))/x/x/x;
	}
	
  return window;
}

inline double WGaussian(double x)
{
  double window;

  window=exp(-x*x/2.0);

  return window;
}

inline double F_Omega(double a)
{
  double omega_a;

  omega_a= Omega/(Omega + a*(1-Omega-OmegaLambda) + a*a*a*OmegaLambda/pow(a,3.0*(1.0+w)));

  return pow(omega_a, 0.6);
}

inline double Jenkins(double sig)
{
  double res,s;

  s = log(1./sig);
  res = 0.315*exp(-pow(fabs(s+0.61),3.8));      

  return res;
}



inline double Sheth_Tormen(double nu)
{
  double res,ca,nu_prime,p;

  ca=0.707;
  p=-0.3;

  nu_prime = sqrt(ca)*nu;

  res = 0.3222*sqrt(2.0/M_PI) *nu_prime*exp(-0.5*nu_prime*nu_prime)/(1.+ 1./pow(nu_prime,0.6));

  return res;
}


inline double Press_Schechter(double nu)
{
  double res;

  res=sqrt(2.0/M_PI) * nu *exp(-0.5*nu*nu);

  return res;
}





inline double delta_c_func()
{ 
  double res;
  //res=1.0 - 0.0052*(1.0-Omega_z) - 0.009*pow(1.0-Omega_z,3.0) - 0.01*pow(1.0-Omega_z, 18.0);
	res = 1.0 - 0.00123*log10(Omega_z);
	return res;
}

inline double efn(double x, double rsphere)
{
  double rk,res;

  rk=exp(x);
  res=rk*var3(rk, rsphere);
  
  return res;
}


inline double var3(double x, double rsphere)
{
  double rk,res,pspec,xk;

  //xk=x/1000.0; // kpc -> Mpc scaling 
  //res=dweight(x, rsphere)*PowerSpec(xk) * (1000.0/Norm);
  res=dweight(x, rsphere)*PowerSpec(x) * (1./Norm);

  return res;
}

inline double var2(double x, double rsphere)
{
  double res,pspec,xk;
  
  //xk=x/1000.0; // kpc -> Mpc scaling 
  //res=weight(x, rsphere)*PowerSpec(xk) * (1000.0/Norm); //renormalize for mf code
  res=weight(x, rsphere)*PowerSpec(x) * (1./Norm);

  return res;
}


inline double evar2(double x, double rsphere)
{
  double rk,res;

  rk = exp(x);
  res = var2(rk, rsphere)*rk;
  
  return res;
}

inline double weight(double x, double rsphere) 
{
  // Tophat filter * k^2
  // there appears to be appreciable differences between C-version and fortran
  double y,res,yinv,yinv2;
  
  y=rsphere*x;
  yinv=1.0/y;
  yinv2=yinv*yinv;

  res=36.0*M_PI*x*x* pow((sin(y)/y-cos(y))*yinv2, 2);
  
  if(sin(y)/y - cos(y) < TINY)
    res=0;

  return res;
}


inline double dweight(double x, double rsphere)  // derivative of weight with rsphere
{
  double y,yinv,yinv2,yinv3,yinv4,res;

  y=rsphere*x;
  yinv=1.0/y;
  yinv2=yinv*yinv;
  yinv3=yinv*yinv2;
  yinv4=yinv*yinv3;


  res=M_PI*x*x*x*72.0*(3.0*cos(y)*yinv3 - 3.0*sin(y)*yinv4+sin(y)*yinv2)
    *(sin(y)*yinv3-cos(y)*yinv2);

  return res;
}


double dndlogm(double logm)
{
  // Number density per interval in log_{10} in mass
  
  double rho_crit=2.7754e11;
  double delta=1.68647;
  double result,rm,rmass,sig,ans,fract;
  double r, res;
  
        
  rm = pow(10,logm);
  rmass = rm/rho_crit/Omega;
  
  sig = sigma_m(rmass, &r);

  ans = dlogdsigma(rmass, r, sig);

  fract = fract_mass(sig);
  if(OPT==4){
    double nu=delta_c_func()*delta/sig;
    double neff=-6*ans-3;

    fract=fract*exp(-0.03/pow(neff+3,2)/pow(nu,0.6));
  }

  res = Rlog_e10*(rho_crit*Omega/rm*fabs(ans)*fract);


  cout<<"rm="<<rm<<" ans="<<ans<<" sig="<<sig<<" fract="<<fract<<endl;


  return res;
}


double sigma_m(double m, double *rsphere_return)
{
  //   Use unit of mass where 1h^{-1}Mpc^3 has mass 1

  double res,rsphere,ling,aend;

  rsphere = pow((3.*m/4./M_PI), 0.33333333333);
  *rsphere_return = rsphere;

  res = sqrt(Const_MF * unnsigma(rsphere));
  
  return res;
}



double dlogdsigma(double mass, double rsphere, double sigma) //dlog(sigma)/dlogM
{
  double rho_crit=2.7754e11;
  double usig,ling,res,aend;

  res=sigdsigdr(rsphere);

  res=Const_MF * res * (mass/sigma/sigma)/4.0/M_PI/rsphere/rsphere;

  return res;
}

double sigdsigdr(double rsphere) //d(sigma^2)/dr
{
  int i;
  double dxk=0.005,xk;
  double sum=0,res;
  
  for(i=0;i<=8000;i++){
    xk=-20.0 + dxk*i;
    sum += efn(xk, rsphere)*dxk;
  }
  
  res = 0.5*sum;
  return res;
}


double  fract_mass(double sig)  //  f(\ln\sigma^{-1})
{
  double delta=1.68647;
  double sqrt_two_over_pi=0.79788456;
  double sqrt_two=1.414213562;
  double nu,nu_prime,s,fract,fdeltac;
  int iselect, ivalid;

  fdeltac=delta_c_func();

  nu = fdeltac*delta/sig;
  
  switch(OPT)
    {
	  case 0:{ // Jenkins et al 2000 fitting formula
      //fract=Jenkins(sig);

      //s = log(1./sig);
      //fract = 0.315*exp(-pow(fabs(s+0.61),3.8));      

      //if (s > -1.2 &&  s < 1.05)
      //	ivalid = 1;
      //else
      //	ivalid=0;
      //break;
      
      //Jenkins et al for FOF0.164 LCDM
      s = log(1./sig);
      fract = 0.301*exp(-pow(fabs(s+0.64),3.88));      

      if (s > -1.2 &&  s < 1.05)
	ivalid = 1;
      else
	ivalid=0;
      break;
		}
	  case 1:{ // Sheth-Tormen
      //fract=Sheth_Tormen(nu);

      nu_prime = sqrt(0.707)*nu;
      fract = 0.3222*sqrt_two_over_pi*nu_prime*exp(-nu_prime*nu_prime/2.)
	*(1.+ 1./pow(nu_prime,0.6));

      ivalid = 1;
      break;
		}
	  case 2:{ // Press-Schechter
      //fract=Press-Schechter(nu);

      fract = sqrt_two_over_pi * nu * exp(-nu*nu/2.);
      ivalid = 1; 
      break;
		}
	  case 3:{ // Yahagi et al;
      
      fract = 0.298*(1.0+pow(0.893*nu/sqrt_two, 1.39))*pow(nu,0.408)
	*exp(-0.893*0.893*nu*nu/2.);
	
      ivalid = 1; 
      break;
		}
	  case 4:{ // Reed et al 2007 fitting formula
      
      s=log(1./sig);
      nu_prime=sqrt(0.707)*nu;
      double G1=exp(-pow((s-0.40)/0.6,2.0)/2);
      double G2=exp(-pow((s-0.75)/0.2,2.0)/2);
      fract= 0.3222*sqrt(2/PI)*nu_prime*exp(-nu*nu*0.746/2.)
	*(1.+1./pow(nu_prime,0.6)+0.6*G1+0.4*G2);
      
      if(s>-1.2 && s<1.05)
	ivalid=1;
      else
	ivalid=0;
      break;
		}
	  case 5:{ // Bhattacharya et ak 2011 fitting formula
		double A,B,p,q;
		A = 0.333/pow(1+calc_z,0.11);
		B = 0.788/pow(1+calc_z,0.01);
		p = 0.807;
		q = 1.795;
		
		fract = A*sqrt(2/PI)*exp(-0.5*B*nu*nu)
				*(1+pow(1/B/nu/nu,p))*pow(sqrt(B)*nu,q);
		ivalid = 1;
		break;
		}	
	}

  return fract;
}




double unnsigma(double rsphere)
{
  int i,j,k;
  double dxk=0.01,xk;
  double sum=0;
  
  for(i=0;i<=4000;i++){
    xk=-20.0 + dxk*i;
    sum += evar2(xk, rsphere)*dxk;
  }
  
  return sum;
}

double RungeKutta(double a_in,double a){
	// u=D/a ,initial condition---> du/dlna=0,u=1
	int i,j;
	double h=(log(a)-log(a_in))/10000;
	double x=log(a_in);
	double u=1;
	double dudlna=0;
	double k0[2],k1[2],k2[2],k3[2];
	
	if(a_in==0){
		printf("you cannot solve calculate linear density growth eq.");
	}if(a == a_in){
		u=1;
	}else{
		for(i=0;i<10000;i++){

		k0[0]=h*dudlna;
		k0[1]=h*(-coeff1(exp(x))*dudlna-coeff2(exp(x))*(u));
			
		k1[0]=h*(dudlna+k0[1]/2);
		k1[1]=h*(-coeff1(exp(x+h/2))*(dudlna+k0[1]/2)-coeff2(exp(x+h/2))*(u+k0[0]/2));
			
		k2[0]=h*(dudlna+k1[1]/2);
		k2[1]=h*(-coeff1(exp(x+h/2))*(dudlna+k1[1]/2)-coeff2(exp(x+h/2))*(u+k1[0]/2));
			
		k3[0]=h*(dudlna+k2[1]);
		k3[1]=h*(-coeff1(exp(x+h))*(dudlna+k2[1])-coeff2(exp(x+h))*(u+k2[0]));
		
		u = u + (k0[0]+2*k1[0]+2*k2[0]+k3[0])/6;
		dudlna = dudlna + (k0[1]+2*k1[1]+2*k2[1]+k3[1])/6;
		x = x+h;
	  }
	}
	
	return a*u;
}


double correction_fac(double logm){
	double rho_crit=2.7754e11;
	double delta=1.68647;
	double result,rm,rmass,sig,ans,fract;
	double r, res;
	double nu,fdeltac;
	
	if(fnl==0){
		res=1;
	}else{
		rm = pow(10,logm);
		rmass = rm/rho_crit/Omega;
  
		sig = sigma_m(rmass, &r);

		ans = dlogdsigma(rmass, r, sig);

		fdeltac= sqrt(0.8)*delta_c_func()*delta;
		nu =fdeltac/sig/Dplus;
	
		double skew=skewness(r)/*1e5*/;
		double dskewdlnsig=dskewdlogr(r)/*1e5*//ans/3;
		double dsigS3dlnsig=-(skew-dskewdlnsig)/pow(sig,3);
		double dS3dlnsig=dsigS3dlnsig/sig - skew/pow(sig,4);
	
	//printf("S3 times sigma^4 =%g\n",S3);
	//printf("d(skewness)/dlogM =%g\n",dSdlnsig*ans);
		
		switch(OPT_correc)
		{
		  case 1:{
				res= 1+(nu*nu*nu-2*nu-1/nu)*skew/pow(sig,3)/6+(nu-1/nu)*dsigS3dlnsig/6;
				break;
			}
		
		  case 2:{
				res = exp(skew*(nu/sig)*(nu/sig)*(nu/sig)/6)
					*fabs(dS3dlnsig*fdeltac/sqrt(1-fdeltac*skew/3/pow(sig,4))/6+sqrt(1-fdeltac*skew/3/pow(sig,4)));
				break;
			}
		
		  default:
			fprintf(stdout,"Not supported\n"); 
			exit(1);
		}
		
	}
	return res;
}


double alpha(double k){
	double res;
	
	res = 2*k*k*tk_eh(k)*C*C/(3*Omega*10000)/Dplus;
	
	return res;
}

double skewness(double r){
	double res,err;
	double xl[3] ={-20,-20,-1};
	double xu[3]={20,20,1};
	const gsl_rng_type *T;
	gsl_rng *random;
	gsl_monte_function F={&skewint,3,&r};
	size_t calls=1000000;
	
	gsl_rng_env_setup();
	
	T=gsl_rng_default;
	random=gsl_rng_alloc(T);
	
	gsl_monte_miser_state *state=gsl_monte_miser_alloc(3);
	gsl_monte_miser_integrate(&F,xl,xu,3,calls,random,state,&res,&err);
	gsl_monte_miser_free(state);
	
	res *= fnl/(2*PI*PI)/(2*PI*PI);
	
	return res;
	
}

double skewint(double *k,size_t dim,void *params){
	double res;
	double r = *(double *)params;
	double k0=exp(k[0]);
	double k1=exp(k[1]);
	//double k0=1/k[0]-1;
	//double k1=1/k[1]-1;
	//double A=Norm*(1.5*Omega*10000/C/C)*(1.5*Omega*10000/C/C);
	double A=0.36*delta_R*(2*PI*PI);
	double q=sqrt(k0*k0+k1*k1+2*k[2]*k0*k1);
	
	res = alpha(k0)*alpha(k1)*alpha(q)*A*pow(k0*HubbleParam/0.002,ns-1.0)/pow(k0,3)*A*pow(k1*HubbleParam/0.002,ns-1.0)/pow(k1,3)*WTopHat(k0*r)*WTopHat(k1*r)*WTopHat(q*r)*(1+2*(k1*k1*k1)/(q*q*q));
	
	return res*(k0)*(k0)*(k1)*(k1)*exp(k[0])*exp(k[1]);
	//return res*(k0)*(k0)*(k1)*(k1)/k[0]/k[0]/k[1]/k[1];
}
	
double dskewdlogr(double r){
	double res,err;
	double xl[3] ={-20,-20,-1};
	double xu[3]={20,20,1};
	const gsl_rng_type *T;
	gsl_rng *random;
	gsl_monte_function F={&skewint2,3,&r};
	size_t calls=1000000;
	
	gsl_rng_env_setup();
	
	T=gsl_rng_default;
	random=gsl_rng_alloc(T);
	
	gsl_monte_miser_state *state=gsl_monte_miser_alloc(3);
	gsl_monte_miser_integrate(&F,xl,xu,3,calls,random,state,&res,&err);
	gsl_monte_miser_free(state);
	
	res *=fnl/(2*PI*PI)/(2*PI*PI);
	
	return res;
	
}

double skewint2(double *k,size_t dim,void *params){
	double res;
	double r = *(double *)params;
	//double A=Norm*(1.5*Omega*10000/C/C)*(1.5*Omega*10000/C/C);
	double A=0.36*delta_R*(2*PI*PI);
	double k0=exp(k[0]);
	double k1=exp(k[1]);
	//double k0=1/k[0]-1;
	//double k1=1/k[1]-1;
	double q=sqrt(k0*k0+k1*k1+2*k[2]*k0*k1);
	
	res = alpha(k0)*alpha(k1)*alpha(q)*A*pow(k0*HubbleParam/0.002,ns-1.0)/pow(k0,3)*A*pow(k1*HubbleParam/0.002,ns-1.0)/pow(k1,3)*(1+2*(k1*k1*k1)/(q*q*q));
	res *= (dWindow(k0*r)*WTopHat(k1*r)*WTopHat(q*r)+WTopHat(k0*r)*dWindow(k1*r)*WTopHat(q*r)+WTopHat(k0*r)*WTopHat(k1*r)*dWindow(q*r));
	
	return res*(k0)*(k0)*(k1)*(k1)*exp(k[0])*exp(k[1]);
	//return res*(k0)*(k0)*(k1)*(k1)/k[0]/k[0]/k[1]/k[1];
}

double dWindow(double x){
	double res;
	res=3*((x*x-3)*sin(x)+3*x*cos(x))/x/x/x;
	return res;
}
