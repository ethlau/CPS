#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "proto.h"
#include "allvars.h"
#include "nrutil.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

using namespace std;

static char TransTable[256];
static double xp[1000],yp[1000],yp2[1000];
static double tab_z[NPOINTS];
static double scale_f[NPOINTS],GF[NPOINTS],GF2[NPOINTS];

static double tab_m[NPOINTS];
static double tab_z_mf[NPOINTS];
static double tab_mf_norm[NPOINTS],err_mf_norm[NPOINTS];
static double **tab_mf,**err_mf;
static double **tab_bh,**err_bh;
static double **tab_cv,**err_cv;
static double **tab_rv,**err_rv;
static double **tab_nuM, **err_nuM;

static int WhichSpectrum, np, WhichWindow, OPT, OPT_fit, OPT_HACK, WhichGF, OPT_offset;
static int OPT_MG, OPT_modified;
static int OPT_survey;

static double bling, scale_Pk;

static double r_tophat,Rsmooth,Delta_c,fdelta_c,Const_MF,st_norm;

static double AA,BB,CC;
static double B1,B2,B3,B4,B5;
static double nu,sigma, Omega_z;

static double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
static double G, Hubble;

static double Norm, InitTime, calc_z, a_exp;
static double Dplus; /* growth factor */

struct get_mvir{
	double delta;
	double rd;
	double redshift;
};

void set_halo_conc_DK15(char *inputPk, double Omega_m, double Omega_b, double w0, double h0, double ns_index){
	
  //input parameters======================================
  WhichSpectrum = 3;  // 1 Eisenstein-Hu, 2 Bond-Efstathiou, 3 CMBFAST table, 4 BBKS
  WhichGF = 2;        // 1 apploximation formula, 2 caluculate linear density growth eq.
  OPT = 6;     // 0 JK  1 ST  2 PS 3 Yahagi 4 Reed et al (astro-ph/0607150v4) 5 Bhattacharya et al (astro-ph/1005.2239) 6 Tinker et al 2008,2010
	
  sprintf(TransTable, "%s", inputPk);
	
  //set cosmological model
  HubbleParam = h0;
  ns = ns_index;
  w  = w0;
  Omega = Omega_m;
  OmegaLambda = 1.-Omega;
  OmegaBaryon = Omega_b;
	
  growth_spline();
	
  int i;
  i=set_units();
  i=initialize_powerspectrum(WhichSpectrum);
	
  stack_table_and_spline();
	
}

void free_halo_conc_DK15(){
	
	free_dmatrix(tab_nuM,1,NPOINTS,1,NPOINTS);
	free_dmatrix(err_nuM,1,NPOINTS,1,NPOINTS);
	free_dmatrix(tab_cv,1,NPOINTS,1,NPOINTS);
	free_dmatrix(err_cv,1,NPOINTS,1,NPOINTS);
	free_dmatrix(tab_mf,1,NPOINTS,1,NPOINTS);
	free_dmatrix(err_mf,1,NPOINTS,1,NPOINTS);
	
}

//==================================================
int initialize_powerspectrum(int Spectrum)
{
	double res;
	int i;
	
	for(i=0;i<1000;i++)
	xp[i]=yp[i]=yp2[i]=0;
	
	if(WhichSpectrum==3){
		//fprintf(stdout,"initialising...\n");
		readCMB_and_do_spline();
	}
	
	a_exp=1/(1+calc_z);
	
	AA=6.4/Gamma; 
	BB=3.0/Gamma; 
	CC=1.7/Gamma;  
	nu=1.13;
	
	
	B1=2.34/Gamma; 
	B2=3.89/Gamma; 
	B3=16.1/Gamma; 
	B4=5.46/Gamma; 
	B5=6.71/Gamma; 
	
	
	Norm = 1.0;
	Sigma8 = sqrt(TopHatSigma2(8.));
	//Norm=Sigma8*Sigma8/res;
	fprintf(stdout,"Sigma8 = %g \n",Sigma8);
	
	Dplus= GrowthFactor(a_exp, 1.0);
	
	return i;
}

//==================================================
int set_units()    /* ... set some units */
{
	UnitLength_in_cm= 3.085678e21; /* 1.0 kpc */
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
double PowerSpec(double kmag)
{
	switch(WhichSpectrum)
	{
	  case 1:
		return PowerSpec_EH(kmag);
		break;
	  case 2:
		return PowerSpec_Efstathiou(kmag);
		break;
	  case 3:
		return PowerSpec_CMBFAST(kmag);
		break;
	  case 4:
		return PowerSpec_BBKS(kmag);
		break;
	  default:
		fprintf(stdout,"Not supported\n");  
	}
}

//==================================================
inline double PowerSpec_Efstathiou(double k)
{
	return Norm*pow(k,ns) / pow(1+pow(AA*k+pow(BB*k,1.5)+CC*CC*k*k,nu),2/nu);
}




//==================================================
inline double PowerSpec_BBKS(double k)
{
	return Norm*pow(k,ns) * pow(log(1.0+B1*k)/(B1*k),2)/ 
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
	return Norm*pow(k,ns)*pow( tk_eh(k), 2);
}



//==================================================
inline double PowerSpec_CMBFAST(double k)
{
	//return Norm*pow(k,ns)*pow(transfunc_cmbfast(k), 2);
	return Norm *transfunc_cmbfast(k);
}


//==================================================
//==================================================
double transfunc_cmbfast(double k)
{
	int i;
	double lk;
	double pow_index;
	
	//k *= 1000.0; /* convert to h/Mpc */
	
	lk=log10(k);
	
	if(lk < xp[0]){
		double dummy = (yp[2]-yp[1])/(xp[2]-xp[1])*(lk-xp[1])+yp[1];
		return pow(10.,dummy);
	} /* usually should not be executed */
	if(lk > xp[np-1]){
		double dummy = (yp[np-2]-yp[np-1])/(xp[np-2]-xp[np-1])*(lk-xp[np-1])+yp[np-1];
		return pow(10.,dummy);
		//return pow(10.,yp[np-1])*pow(k/pow(10.,xp[np-1]),-2.);
		//return pow(10.,yp[np-1]);
	}
	
	splint(xp-1, yp-1, yp2-1, np, lk, &pow_index);
	
	return pow(10.0,pow_index);
}

//==================================================
void readCMB_and_do_spline()
{
	int i,iline;
	double yp1,ypn;
	char tfunc_table[100];
	int  errorFlag=0;
	FILE *fd;
	
	sprintf(tfunc_table,"%s",TransTable);
	
	fprintf(stdout,"Reading %s .\n",tfunc_table);
	iline=0;
	double dummy1,dummy2,dummy3,dummy4,dummy5;
	if(fd=fopen(tfunc_table,"r")){
		while(!feof(fd)){
			fscanf(fd,"%lf %lf\n",&xp[iline],&yp[iline]);
			//fprintf(stdout,"%g %g \n",xp[iline],yp[iline]);
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
	
	//for(i=0;i<iline;i++)
	//fprintf(stdout,"%g %g \n",xp[i],yp[i]);
	
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
	double kr,kr3,kr2,wf,x;
	
	kr=r_tophat*k;
	kr2=kr*kr;
	kr3=kr2*kr;
	
	if(kr<1e-8) return 0;
	
	wf=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
	x=4*PI*k*k*wf*wf*PowerSpec(k)/pow(2*PI,3.);
	
	return x;
}

//==================================================
double GrowthFactor(double astart, double aend)
{
	return growth(aend)/growth(astart);
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
		return growth_for_any_w(a);
		break;
	  default:
		fprintf(stdout,"Not supported\n"); 
	}
}

inline double growth_int(double a)
{
	return pow(a / (Omega + (1-Omega-OmegaLambda)*a + OmegaLambda*a*a*a), 1.5);
}


inline double var2(double x, double rsphere)
{
	double res,pspec,xk;
	res=weight(x, rsphere)*PowerSpec(x)*bling*bling;
	
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
	
	res=36.0*M_PI*x*x* pow((sin(y)/y-cos(y))*yinv2, 2)/pow(2*PI,3);
	
	if(sin(y)/y - cos(y) < TINY)
	res=0;
	
	return res;
}

double sigma_m(double m, double *rsphere_return)
{
	//   Use unit of mass where 1h^{-1}Mpc^3 has mass 1
	
	double res,rsphere,ling,aend;
	
	rsphere = pow((3.*m/4./M_PI), 0.33333333333);
	*rsphere_return = rsphere;
	
	res = sqrt(unnsigma(rsphere));
	
	return res;
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

double delta_v(double z){
	double Omz = (Omega*pow(1.0+z,3.0))/(Omega*pow(1.0+z,3.0)+OmegaLambda);
	return 18.0*M_PI*M_PI+82.*(Omz-1)-39.0*(Omz-1)*(Omz-1);
}

double r_vir(double z,double M){
	double rho_crit,logr;
	rho_crit=2.7754e11*(Omega*pow(1.0+z,3.0) + OmegaLambda);
	logr=(log10(3*M)-log10(4.0*PI*rho_crit*delta_v(z)))/3.0;
	return pow(10.0,logr);
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

void growth_spline(){
	int i;
	double yp1,ypn;
	double da = (1.0-(1.0/(1.0+1088.2)))/NPOINTS;
	
	for(i=0;i<NPOINTS;i++){
		scale_f[i] = (double)(i+1)*da + 1.0/(1.0+1088.2);
		GF[i] = RungeKutta(1.0/(1.0+1088.2),scale_f[i]);
	}
	
	for(i=0;i<NPOINTS;i++){
		scale_f[i] = log10(scale_f[i]);
		GF[i] = log10(GF[i]);
	}
	
	yp1 = 1.e31;
	ypn = 1.e31;
	
	spline(scale_f-1, GF-1, NPOINTS, yp1, ypn, GF2-1);
	
}

double growth_for_any_w(double a){
	double la;
	double pow_index;
	
	la=log10(a);
	
	if(la < scale_f[0])
	return a; /* usually should not be executed */
	if(la > scale_f[NPOINTS-1]){
		double dummy = (GF[NPOINTS-2]-GF[NPOINTS-1])/(GF[NPOINTS-2]-GF[NPOINTS-1])*(la-scale_f[NPOINTS-1])+GF[NPOINTS-1];
		return pow(10.,dummy);
		//return pow(10.,yp[np-1]);
	}
	
	splint(scale_f-1, GF-1, GF2-1, NPOINTS, la, &pow_index);
	return pow(10.0,pow_index);
}

void stack_table_and_spline(){
  int i,j,k;
  double yp1,ypn;
  double dlogz = (log10(4.0)-(-3.))/(NPOINTS-1);
  double dlogm = (16.5-10.)/(NPOINTS-1);
	
  yp1 = 1.e31;
  ypn = 1.e31;
	
  for(i=0;i<NPOINTS;i++){
    tab_z[i] = i*dlogz + (-3.);
  }
	
  tab_nuM=dmatrix(1,NPOINTS,1,NPOINTS);
  err_nuM=dmatrix(1,NPOINTS,1,NPOINTS);
  tab_cv=dmatrix(1,NPOINTS,1,NPOINTS);
  err_cv=dmatrix(1,NPOINTS,1,NPOINTS);
	
  for(i=0;i<NPOINTS;i++){
    a_exp=1.0/(1.0+pow(10,tab_z[i])); // expansion parameter at which MF is computed
    calc_z=1.0/a_exp -1.0;
		
    Delta_c=1.68647;
		
    Omega_z = Omega*pow(1.0+calc_z,3.0)/(Omega*pow(1.0+calc_z,3.0) + OmegaLambda);
    bling = growth(a_exp)/growth(1.);
		
    for(j=0;j<NPOINTS;j++){
      tab_m[j] = j*dlogm + 10.;
      tab_nuM[i+1][j+1] = nu_M(calc_z, pow(10., tab_m[j]));
    }
  }
	
  for(i=1;i<=NPOINTS;i++){
    for(j=1;j<=NPOINTS;j++){
      tab_nuM[i][j] = log10(tab_nuM[i][j]);
    }
  }
  splie2(tab_z-1, tab_m-1,tab_nuM,NPOINTS,NPOINTS,err_nuM);
  printf("set spline nu(z, M)...\n");
	
  for(i=0;i<NPOINTS;i++){
    for(j=0;j<NPOINTS;j++){
			
      //convert Mvir to M200c
      double mvir  = pow(10., tab_m[j]);
      double m200c = M_vir_to_M_200c(pow(10,tab_z[i]), mvir);
			
      double rvir = r_vir(pow(10,tab_z[i]), mvir);
      double r200c = r_delta(pow(10,tab_z[i]), m200c, -200.);
			
      tab_cv[i+1][j+1] = c_200c_DK15(pow(10,tab_z[i]), m200c) * rvir/r200c;
      //cout << pow(10,tab_z[i]) << " " << pow(10., tab_m[j]) << " " << m200c << " " << tab_cv[i+1][j+1] << endl;
      tab_cv[i+1][j+1] = log10(tab_cv[i+1][j+1]);
    }
  }
	
  splie2(tab_z-1, tab_m-1,tab_cv,NPOINTS,NPOINTS,err_cv);
  printf("set spline c_vir(z, M)...\n");

  tab_mf=dmatrix(1,NPOINTS,1,NPOINTS);
  err_mf=dmatrix(1,NPOINTS,1,NPOINTS);
	
  for(i=0;i<NPOINTS;i++){

    a_exp=1.0/(1.0+pow(10,tab_z[i])); // expansion parameter at which MF is computed
    calc_z=1.0/a_exp -1.0;
		
    Delta_c=1.68647;
		
    Omega_z = Omega*pow(1.0+calc_z,3.0)/(Omega*pow(1.0+calc_z,3.0) + OmegaLambda);
    bling = growth(a_exp)/growth(1.);
		
    for(j=0;j<NPOINTS;j++){
      tab_mf[i+1][j+1] = dndlogm(tab_m[j]);
      if(tab_mf[i+1][j+1] < 1e-100){
	printf("z=%e M=%e dn/dlogM= %e\n",pow(10,tab_z[i]),pow(10,tab_m[j]),tab_mf[i+1][j+1]);
      }
    }
  }
  //exit(1);
	
  for(i=1;i<=NPOINTS;i++){
    for(j=1;j<=NPOINTS;j++){
      tab_mf[i][j] = log10(tab_mf[i][j]);
    }
  }		
  splie2(tab_z-1, tab_m-1,tab_mf,NPOINTS,NPOINTS,err_mf);
  printf("set spline dn/dlogm...\n");
	
}

void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a){
	int j;
	for(j=1;j<=m;j++){
		spline(x2a,ya[j],n,1.e31,1.e31,y2a[j]);
	}
}
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y){
	int j;
	double *ytmp,*yytmp;
	
	ytmp=dvector(1,n);
	yytmp=dvector(1,n);
	
	for(j=1;j<=m;j++){
		splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	}
	spline(x1a,yytmp,m,1.e31,1.e31,ytmp);
	splint(x1a,yytmp,ytmp,m,x1,y);
	free_dvector(yytmp,1,n);
	free_dvector(ytmp,1,n);
}

double dlnP_dlnk(double lnk){
	double dlnk = 0.01;
	double k1 = exp(lnk-dlnk);
	double k2 = exp(lnk+dlnk);
	double P1 = PowerSpec_EH(k1);
	double P2 = PowerSpec_EH(k2);
	return log(P2/P1)/dlnk/2;
}

double nu_M(double z, double M){
	
	double rho_crit=2.7754e11;
	double r, rmass = M/rho_crit/Omega;
	double sig = sigma_m(rmass, &r);
	
	return 1.686/sig;
}

double nu_M_fast(double z, double M){
	double logm = log10(M);
	double lz = log10(z);
	double pow_index;
	
	splin2(tab_z-1, tab_m-1, tab_nuM, err_nuM, NPOINTS, NPOINTS, lz, logm, &pow_index);
	//printf("OK.spline dn/dlogm\n");
	if( pow_index != pow_index ){
		printf("fail to spline nu at m=%e z=%e\n",pow(10.,logm),z);
		exit(1);
	}
	return pow(10.0,pow_index);
}

double c_200c_DK15(double z, double M){ //https://arxiv.org/pdf/1407.4730.pdf
	double nu = nu_M_fast(z, M);
	double rho_crit=2.7754e11;
	double R = pow(3*M/rho_crit/Omega/(4*M_PI), 0.3333);
	
	//paramaters
	double kap = 0.69;
	double phi0 = 7.14;
	double phi1 = 1.60;
	double eta0 = 4.10;
	double eta1 = 0.75;
	double alpha = 1.40;
	double beta = 0.67;
	
	double n = dlnP_dlnk(log(2*M_PI/R*kap));
	
	double cmin = phi0 + phi1 *n;
	double numin = eta0 + eta1 *n;
	
	double res = cmin/2 * (pow(nu/numin, -alpha) + pow(nu/numin, beta));
	//cout << n << " " << res << endl;
	return res;
}

double c_vir_DK15_fast(double z, double M){
	double logm = log10(M);
	double lz = log10(z);
	double pow_index;
	
	splin2(tab_z-1, tab_m-1, tab_cv, err_cv, NPOINTS, NPOINTS, lz, logm, &pow_index);
	//printf("OK.spline dn/dlogm\n");
	if( pow_index != pow_index ){
		printf("fail to spline cvir at m=%e z=%e\n",pow(10.,logm),z);
		exit(1);
	}
	return pow(10.0,pow_index);
}

double get_eq_for_M_200c(double x, void *params){
	struct get_mvir *p = (struct get_mvir *)params;
	
	double delta = p->delta;
	double rd = p->rd; 
	double z = p->redshift;
	
	double rhom =2.7754e11*Omega*pow(1+z,3);
	double rhoc =2.7754e11*(Omega*pow(1+z,3)+OmegaLambda);
	
	//solve 3 rho_s (z, m200c) rs(z, m200c)^3 f(rd/rs(z, m200c)) = rho_c * delta_vir * r_vir^3 where f(x) = ln(1+x)-x/(1+x)	
	double m200c = pow(10., x);
	double r200c = pow(10., (log10(3*m200c)-log10(4.0*PI*rhoc*200.))/3.0);
	double c200c = c_200c_DK15(z, m200c);
	double rs = r200c/c200c;
	double rhos = m200c/(log(1.+c200c)-c200c/(1.+c200c))/(4*PI*rs*rs*rs);
	double cd  = rd/rs;
	
	double rhs = 3 * rhos * rs* rs* rs* (log(1+cd) -cd/(1+cd));
	double lfs = rhoc * delta * rd * rd * rd;
	return rhs-lfs;
}

double M_vir_to_M_200c(double z, double mvir){
	double rd = r_vir(z, mvir);
	
	//get mvir
	int status; 
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 1.0;
	double x_lo = 5.0, x_hi = 20.0;
	gsl_function F;
	struct get_mvir params={delta_v(z), rd, z};
	
	//printf("%e %e\n",get_eq_for_M_vir(x_lo, &params), get_eq_for_M_vir(x_hi, &params));
	
	F.function = &get_eq_for_M_200c;
	F.params = &params; 
	T = gsl_root_fsolver_brent; 
	s = gsl_root_fsolver_alloc(T);
	
	gsl_root_fsolver_set(s, &F, x_lo, x_hi); 
	
	do {
		iter++;
		status = gsl_root_fsolver_iterate(s); 
		r = gsl_root_fsolver_root(s); 
		x_lo = gsl_root_fsolver_x_lower(s); 
		x_hi = gsl_root_fsolver_x_upper(s); 
		status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-7); 
		//if (status == GSL_SUCCESS) printf ("Converged:\n");
		//printf("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r,  x_hi-x_lo);
		
	} while (status == GSL_CONTINUE && iter < max_iter);
	
	gsl_root_fsolver_free(s);
	
	double m200c = pow(10., r);
	return m200c;
}

// halo mass function

inline double delta_c_func()
{ 
	double res;
	res = 3.0*pow(12*PI,0.666666666)*(1.0+0.00123*log10(Omega_z))/20.0;
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
	
	res=dweight(x, rsphere)*PowerSpec(x) * bling * bling;
	
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
	
	return res/pow(2*PI, 3);
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
  if(OPT == 4){
    double nu=delta_c_func()/sig;
    double neff=-6*ans-3;
		
    fract=fract*exp(-0.03/pow(neff+3,2)/pow(nu,0.6));
  }
	
  res = Rlog_e10*(rho_crit*Omega/rm*fabs(ans)*fract);
		
  //cout<<"rm="<<rm<<" ans="<<ans<<" sig="<<sig<<" fract="<<fract<<endl;
		
  return res;
}

double dlogdsigma(double mass, double rsphere, double sigma) //dlog(sigma)/dlogM
{
	double rho_crit=2.7754e11;
	double usig,ling,res,aend;
	
	res=sigdsigdr(rsphere);
	
	res= res * (mass/sigma/sigma)/4.0/M_PI/rsphere/rsphere;
	
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
  double nu_prime,s,fract,fdeltac;
  int iselect, ivalid;
	
  fdeltac=delta_c_func();
	
  double Nu = fdeltac/sig;
	
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
			
      nu_prime = sqrt(0.707)*Nu;
      fract = 0.3222*sqrt_two_over_pi*nu_prime*exp(-nu_prime*nu_prime/2.)
	*(1.+ 1./pow(nu_prime,0.6));
			
      ivalid = 1;
      break;
    }
    case 2:{ // Press-Schechter
      //fract=Press-Schechter(nu);
			
      fract = sqrt_two_over_pi * Nu * exp(-Nu*Nu/2.);
      ivalid = 1; 
      break;
    }
    case 3:{ // Yahagi et al;
			
      fract = 0.298*(1.0+pow(0.893*Nu/sqrt_two, 1.39))*pow(Nu,0.408)
	*exp(-0.893*0.893*Nu*Nu/2.);
			
      ivalid = 1; 
      break;
    }
    case 4:{ // Reed et al 2007 fitting formula
			
      s=log(1./sig);
      nu_prime=sqrt(0.707)*Nu;
      double G1=exp(-pow((s-0.40)/0.6,2.0)/2);
      double G2=exp(-pow((s-0.75)/0.2,2.0)/2);
      fract= 0.3222*sqrt(2/PI)*nu_prime*exp(-Nu*Nu*0.746/2.)
	*(1.+1./pow(nu_prime,0.6)+0.6*G1+0.4*G2);
			
      if(s>-1.2 && s<1.05)
	ivalid=1;
      else
	ivalid=0;
      break;
    }
    case 5:{ // Bhattacharya et al 2011 fitting formula
      double A,B,p,q;
      A = 0.333/pow(1+calc_z,0.11);
      B = 0.788/pow(1+calc_z,0.01);
      p = 0.807;
      q = 1.795;
			
      fract = A*sqrt(2/PI)*exp(-0.5*B*Nu*Nu)
	*(1+pow(1/B/Nu/Nu,p))*pow(sqrt(B)*Nu,q);
      ivalid = 1;
      break;
    }
    case 6:{ 
      // Tinker et al 2008
      double A,a,b,c;
      double Delta=200.0;
      double alpha = pow(10.,-pow(0.75/(log10(Delta/75.)),1.2));
      A = 0.186*pow(1.+calc_z, -0.14);
      a = 1.47*pow(1.+calc_z,-0.06);
      b = 2.57*pow(1.+calc_z, -alpha);
      c = 1.19;
      fract = A*(1+pow(sig/b,-a))*exp(-c/sig/sig);
			
      // Tinker et al 2010 (http://arxiv.org/pdf/1001.3162v2.pdf)
      /*
	Nu = 1.69/sig;
	double alp, bet, gam, phi, eta;
	alp = 0.368;
	bet = 0.589;
	gam = 0.864;
	phi = -0.729;
	eta = -0.243;
			
	bet = bet * pow(1.+calc_z, 0.20);
	gam = gam * pow(1.+calc_z, -0.01);
	phi = phi * pow(1.+calc_z, -0.08);
	eta = eta * pow(1.+calc_z, 0.27);
			
	fract = alp * Nu * (1+pow(bet*Nu, -2*phi))*pow(Nu, 2*eta)*exp(-0.5*gam*Nu*Nu);
      */
      ivalid = 1;
      break;
    }
    }
	
  return fract;
}

double dndlogm_fast(double logm, double z){
  double lz = log10(z);
  double pow_index;
	
  splin2(tab_z-1, tab_m-1, tab_mf, err_mf, NPOINTS, NPOINTS, lz, logm, &pow_index);
  //printf("OK.spline dn/dlogm\n");
  if( pow_index !=pow_index ){
    printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
    exit(1);
  }
  return pow(10.0,pow_index);
}

// Mvir <-> Mdelta for Tinker mass function / halo bias
double r_delta(double z, double Mass, double delta){
  double rho_crit,logr;
	
  if(delta > 0){
    //This difinition is the same as that in Tinker et al
    double rhom=2.7754e11*Omega*pow(1+z,3);
    logr=(log10(3*Mass)-log10(4.0*PI*rhom*fabs(delta)))/3.0;
    return pow(10.0,logr);
  }else{
    //This difinition is the same as that in Prada et al
    rho_crit = 2.7754e11*(Omega*pow(1+z,3)+OmegaLambda);
    logr=(log10(3*Mass)-log10(4.0*PI*rho_crit*fabs(delta)))/3.0;
    return pow(10.0,logr);
  }
}

double get_eq_for_M_delta(double x, void *params){
	struct get_mvir *p = (struct get_mvir *)params;
	
	double delta = p->delta;
	double mvir = p->rd; //this is mvir, here
	double z = p->redshift;
	
	//solve delta/delta_vir f(1/c) = f(1/c rvir/rdelta) where f(x)=x^3 [ln(1+1/x)-1/(1+x)]
	double cv = c_vir_DK15_fast(z, mvir);
	double rs = r_vir(z, mvir)/cv;
	double cd = pow(10.,x);
	
	//printf("%e %e %e %e\n",delta, delta_v(z), cd, cv);
	
	double rhs = delta*Omega*pow(1.+z, 3)/delta_v(z)/(Omega*pow(1+z,3)+OmegaLambda)*(1./cv)*(1./cv)*(1./cv)*(log(1+cv)-cv/(1+cv));
	double lfs = (1./cd)*(1./cd)*(1./cd)*(log(1+cd)-cd/(1+cd));
	return rhs-lfs;
}

double M_vir_to_M_delta(double z, double mvir, double delta){
  //get mvir
  int status; 
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 12.;
  double x_lo = -5.0, x_hi = +5.;
  gsl_function F;
  struct get_mvir params={delta, mvir, z};
	
  //printf("%e %e\n",get_eq_for_M_vir(x_lo, &params), get_eq_for_M_vir(x_hi, &params));
	
  F.function = &get_eq_for_M_delta; 
  F.params = &params; 
  T = gsl_root_fsolver_brent; 
  s = gsl_root_fsolver_alloc(T);
	
  gsl_root_fsolver_set(s, &F, x_lo, x_hi); 
	
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s); 
    r = gsl_root_fsolver_root(s); 
    x_lo = gsl_root_fsolver_x_lower(s); 
    x_hi = gsl_root_fsolver_x_upper(s); 
    status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-7); 
    //if (status == GSL_SUCCESS) printf ("Converged:\n");
    //printf("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r,  x_hi-x_lo);
		
  } while (status == GSL_CONTINUE && iter < max_iter);
	
  gsl_root_fsolver_free(s);
	
  double cd = pow(10.,r);
  double rs = r_vir(z, mvir)/c_vir_DK15_fast(z, mvir);
  double rd = cd*rs;
  double rho_crit=2.7754e11*(Omega*pow(1+z,3)/*+OmegaLambda*/);
  double mdelta = 4*PI/3*delta*rho_crit*rd*rd*rd;
	
  return mdelta;
}
