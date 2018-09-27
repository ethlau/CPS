#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <vector>

#include "cosmo.h"
#include "cluster.h"
#include "gas_model.h"

#include "allvars.h"
#include "nrutil.h"
#include "proto.h"

#include "cfortran.h"

#define Nx 100

static double xp[Nx], yp[Nx], yp2[Nx];
static int Nspline;

using namespace std;

vector<double> calc_Shaw_pressure_profile(cosmo cosm_model, float z, float Mvir, vector<float> x);

struct Shaw_param{
  double alpha0; // fiducial : 0.18
  double n_nt;   // fiducial : 0.80
  double beta;   // fiducial : 0.50
  double eps_f;  // fiducial : 3.97e-6
  double eps_DM; // fiducial : 0.00
  double f_star; // fiducial : 0.026
  double S_star; // fiducial : 0.12
  double A_C;    // fiducial : 1.00
};
static struct Shaw_param S;

vector<double> calc_Flender_pressure_profile(cosmo cosm_model, float z, float Mvir, vector<float> x);

struct Flender_param{
  double alpha0; // fiducial : 0.18
  double n_nt;   // fiducial : 0.80
  double beta;   // fiducial : 0.50
  double eps_f;  // fiducial : 3.97e-6
  double eps_DM; // fiducial : 0.00
  double f_star; // fiducial : 0.026
  double S_star; // fiducial : 0.12
  double A_C;    // fiducial : 1.00
  double gamma_mod0; // fiducial : 0.10
  double gamma_mod_zslope; // fiducial : 1.72
  double x_break; // fiducial : 0.195
  double x_smooth; // fiducial : 0.01
  double n_nt_mod;  // fiducial : 0.80
};
static struct Flender_param F;

void free_FFTdata();
void FFT_density_profile(double *output, double *bin, int nbin);

extern "C"
{
  void fhti_(int * n, double * mu, double * q, double * dlnr,
	     double * kr, int * kropt, double * wsave, int * ok);
  void fht_(int * n, double * a , int * dir, double * wsave);
  void fhtq_(int * n, double * a , int * dir, double * wsave);
}

//FFTlog paramaters
static double logrmin;
static double logrmax;
static int N;

static double q;
static double kr;
static int kropt;
//kropt = 0 to use input kr as is;                                                                                                     
//        1 to change kr to nearest low-ringing kr, quietly;                                                                           
//        2 to change kr to nearest low-ringing kr, verbosely;                                                                         
//        3 for option to change kr interactively.
static int dir;
static double *r, *a, *k, *wsave;

void set_FFTlog_param(){
	
  logrmin=log10(1e-8);
  logrmax=log10(1e+8);
  N=512;
	
  kr=1.;
  kropt=1;
	
  dir=1;
	
  r= new double [N];
  a= new double [N];
  k= new double [N];
  wsave = new double [5*N];
	
}

int main(int argc, char **argv){
	
  char outname[250];
	
  if(argc!=2){
    fprintf(stderr,"usage:\n > save_yl input\n");
    return 1;
  }
	
  FILE *fin;
  fin = fopen(argv[1], "r");
  if(fin == NULL){
    fprintf(stderr, "you can not find %s\n", argv[1]);
    exit(1);
  }
	
  /* Cosmological parameters */
  float H0, Omega_M, Omega_b, wt, Omega_k, ns;
  fscanf(fin, "%f %f %f %f %f %f\n", &H0, &Omega_M, &Omega_b, &wt, &Omega_k, &ns);
	
  double X = 0.76;// primordial hydrogen fraction
  double factor = (2.0*X+2.0)/(5.0*X+3.0);//conversion factor from gas pressure to electron pressure
	
  cout << H0 << " " << Omega_M << " " << Omega_b << " " << wt << " " << Omega_k << " " << ns << endl;
	
  /* cosmology class */
  cosmo cosm_model(H0, Omega_M, Omega_b, Omega_k, wt);
	
  char inputPk[256];
  fscanf(fin, "%s", inputPk);
	
  set_halo_conc_DK15(inputPk, Omega_M, Omega_b, wt, H0/100., ns);
	
  /* Shaw et al paramaters */
  //fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &S.alpha0, &S.n_nt, &S.beta, &S.eps_f, &S.eps_DM, &S.f_star, &S.S_star, &S.A_C);
	
  /* Flender et al paramaters */
  fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
	 &F.alpha0, &F.n_nt, &F.beta, &F.eps_f, &F.eps_DM, &F.f_star, &F.S_star, &F.A_C, 
	 &F.gamma_mod0, &F.gamma_mod_zslope, &F.x_break, &F.x_smooth, &F.n_nt_mod);

	
  int nzbin;
  float zmin, zmax;
  fscanf(fin, "%d %f %f\n", &nzbin, &zmin, &zmax);
  cout << nzbin << " " << zmin << " " << zmax << endl;
	
  int nmbin;
  float logMvir_min, logMvir_max;
  fscanf(fin, "%d %f %f\n", &nmbin, &logMvir_min, &logMvir_max);
  cout << nmbin << " " << logMvir_min << " " << logMvir_max << endl;
	
  float z, Mvir, x;
  vector<float> zlist, Mlist, xlist;
  /* set up for redshift, virial mass, radius bin */
  //redshift
  for(int i=0;i<nzbin;i++){
    z = (log(zmax)-log(zmin))/(nzbin-1)*(float)(1.*i) + log(zmin);
    z = exp(z);
    zlist.push_back(z);
  }
	
  //virial mass (Brayn & Norman) in [Msun], not in [Msun/h]
  for(int i=0;i<nmbin;i++){
    Mvir = (logMvir_max-logMvir_min)/(nmbin-1)*(float)(1.*i) + logMvir_min;
    Mvir = pow(10.0, Mvir);
    Mlist.push_back(Mvir);
  }
	
  //radius in unit of R500, x = r/R500
  float xmin = 0.01;
  float xmax = 100.0;
  for(int i=0;i<Nx;i++){
    x = (log10(xmax)-log10(xmin))/Nx*(float)(1.*i+0.5) + log10(xmin);
    x = pow(10.0, x);
    xlist.push_back(x);
  }
	
  //set maximum in integral
  /*
    double xout;
    fscanf(fin, "%lf\n", &xout);
    cout << xout << endl;
  */
	
  //set the output file name
  fscanf(fin, "%s\n", outname);
  cout << outname << endl;
  fclose(fin);

  set_FFTlog_param();
	
  FILE *fp;	
  fp = fopen(outname, "wb");
  if(fp == NULL){
    printf("you can not make %s.\n", outname);
    exit(1);
  }
	
  fwrite(&nzbin, sizeof(int), 1, fp);
  fwrite(&zmin, sizeof(float), 1, fp);
  fwrite(&zmax, sizeof(float), 1, fp);
	
  fwrite(&nmbin, sizeof(int), 1, fp);
  fwrite(&logMvir_min, sizeof(float), 1, fp);
  fwrite(&logMvir_max, sizeof(float), 1, fp);
	
  double dlog_ell = (3.-(-5.))/(Nx-1);
  double tab_l_ls[Nx];
  double bin[Nx];
  double tab_yl_int[Nx];
	
  for(int i=0;i<Nx;i++){
    tab_l_ls[i] = -5.0+(double)(i)*dlog_ell;
    bin[i] = pow(10., tab_l_ls[i]);
  }
		
  int NPOINTS = Nx;
  fwrite(&NPOINTS, sizeof(int), 1, fp);
  fwrite(tab_l_ls, sizeof(double), NPOINTS, fp);
	
  for(int i=0;i<nzbin;i++){
    fprintf(stdout, ".");fflush(stdout);
    for(int j=0;j<nmbin;j++){
			
      vector<double> pressure;
      pressure = calc_Flender_pressure_profile(cosm_model, zlist[i], Mlist[j], xlist);
      //cout << "#x electron pressure[keV/cm^3]" << endl;
      //for(int k=0;k<100;k++) cout << xlist[k] << " " << factor*pressure[k] << endl;
			
      double yp1 = 1.e31;
      double ypn = 1.e31;
			
      Nspline=-1;
      double xout;
      for(int k=0;k<Nx;k++){
	xp[k] = log10(xlist[k]);
	if(factor*pressure[k] < 1e-30){
	  //fprintf(stderr, "found negative pressure at x=%e\n", xlist[k]);
	  //cout << "#z=" << zlist[i] << ",Mvir=" << Mlist[j] << " [Msun]" << endl;
	  //exit(1);
	  xout = pow(10., xp[k-1]);
	  Nspline = k-1;
	  goto NEXT;
	}
	yp[k] = log10(factor*pressure[k]);
      }
    NEXT:;
      if(Nspline < 0){
	fprintf(stderr, "do not find maximum of x = r/r500\n");
	exit(1);
      }
      spline(xp-1, yp-1, Nspline, yp1, ypn, yp2-1);
      FFT_density_profile(tab_yl_int, bin, Nx);	      
      
      /*
      for(int k=0;k<Nx;k++){
	tab_yl_int[k] = y_l_integral(pow(10., tab_l_ls[k]), xout);
	if(isnan(tab_yl_int[k])){
	  fprintf(stderr, "found NaN!\n");
	  cout << "#z=" << zlist[i] << ",Mvir=" << Mlist[j] << " [Msun]" << endl;
	  exit(1);
	}
      }
      */

      fwrite(tab_yl_int, sizeof(double), NPOINTS, fp);
    }
  }
  fclose(fp);
  fprintf(stdout, "DONE!\n");
	
  free_halo_conc_DK15();
  free_FFTdata();

  return 0;
}

vector<double> calc_Shaw_pressure_profile(cosmo cosm_model, float z, float Mvir, vector<float> x){
  /* Shaw model parameters. no need to change so far*/
  float conc_norm = S.A_C;
  float conc_mass_norm = 1.0;
  float ad_index = 5.0; // Gamma = 1+1./ad_index in arXiv:1706.08972
  /*
    float delta_rel = 0.18, delta_rel_n = 0.8, delta_rel_zslope = 0.5; // delta_rel = alpha_0, delta_rel_n  = n_nt, delta_rel_zslope =  beta in Shaw et al 2010
    float eps_fb = 3.97e-6; // epsilon_f in arXiv:1706.08972
    float eps_dm = 0.0; // epsilon_DM in arXiv:1706.08972
    float fs_0 = 0.026; // f_star in arXiv:1706.08972
    float fs_alpha = 0.12; // S_star in arXiv:1706.08972
  */
  float delta_rel = S.alpha0, delta_rel_n = S.n_nt, delta_rel_zslope = S.beta;
  float eps_fb = S.eps_f;
  float eps_dm = S.eps_DM;
  float fs_0 = S.f_star;
  float fs_alpha = S.S_star;
	
  int pturbrad = 2;
  bool verbose = false;
  float overden_id = -1.0;
  int relation = 3;
  float rcutoff = 2.0;
  float Omega_M = cosm_model.get_Omega_M();
  float Omega_b = cosm_model.get_Omega_b();
  float h = cosm_model.get_H0()/100;
  float cosmic_t, cosmic_t0, M500, R500;
  double r, pre;
  vector<double> pressure;
	
  cosmic_t = cosm_model.cosmic_time(z);
  cosmic_t0 = cosm_model.cosmic_time(0.0);
  //cluster nfwclus(Mvir, z, overden_id, relation, cosm_model);
  //nfwclus.concentration(conc_norm, conc_mass_norm); // set halo concentration using M-c relation of Duffy et al (08)
  //M500 = nfwclus.get_mass_overden(500.0);// Msun
  //R500 = nfwclus.get_rad_overden(500.0);// (physical) Mpc
	
  cluster nfwclus(Mvir, z, overden_id, relation, cosm_model);
  float cvir = conc_norm * c_vir_DK15_fast(z, Mvir*h);
	
  nfwclus.set_conc(cvir);
  M500 = nfwclus.get_mass_overden(500.0);// Msun
  R500 = nfwclus.get_rad_overden(500.0);// (physical) Mpc
	
  gas_model icm_mod(delta_rel, ad_index, eps_fb, eps_dm, fs_0, fs_alpha, pturbrad, delta_rel_zslope, delta_rel_n);
	
  icm_mod.calc_fs(M500, Omega_b/Omega_M, cosmic_t0, cosmic_t);
  icm_mod.evolve_pturb_norm(z, rcutoff);
  icm_mod.set_nfw_params(Mvir, nfwclus.get_radius(), nfwclus.get_conc(), nfwclus.get_rhoi(), R500);
  icm_mod.set_mgas_init(Omega_b/Omega_M);
  icm_mod.findxs();
	
  icm_mod.solve_gas_model(verbose, 1e-5);
	
  for(int xi=0;xi<x.size();xi++){
    r = (double) x[xi]*R500;
    pre = icm_mod.calc_gas_pressure(r, R500); // keV/cm^3
    if(pre < 0.0) pre = 0.0;
    pressure.push_back(pre);
  }
	
  return pressure;
}

vector<double> calc_Flender_pressure_profile(cosmo cosm_model, float z, float Mvir, vector<float> x){
	
  float conc_norm = F.A_C;
  float conc_mass_norm = 1.0;
  float ad_index = 5.0; // Gamma = 1+1./ad_index in arXiv:1706.08972
  
  /*
    float delta_rel = 0.18, delta_rel_n = 0.8, delta_rel_zslope = 0.5; // delta_rel = alpha_0, delta_rel_n  = n_nt, delta_rel_zslope =  beta in Shaw et al 2010
  
    float eps_fb = 3.97e-6; // epsilon_f in arXiv:1706.08972
    float eps_dm = 0.0; // epsilon_DM in arXiv:1706.08972
    float fs_0 = 0.026; // f_star in arXiv:1706.08972
    float fs_alpha = 0.12; // S_star in arXiv:1706.08972
  */
  float delta_rel = F.alpha0, delta_rel_n = F.n_nt, delta_rel_zslope = F.beta;
  float eps_fb = F.eps_f;
  float eps_dm = F.eps_DM;
  float fs_0 = F.f_star;
  float fs_alpha = F.S_star;

  float gamma_mod0 = F.gamma_mod0;
  float gamma_mod_zslope = F.gamma_mod_zslope;
  float x_break = F.x_break;
  float x_smooth = F.x_smooth;

  int pturbrad = 2;
  bool verbose = false;
  float Rvir, M500, R500, Rscale, conc, cosmic_t, cosmic_t0;
  float Omega_M = cosm_model.get_Omega_M();
  float Omega_b = cosm_model.get_Omega_b();
  float h =cosm_model.get_H0()/100.0;
  float E;
  // set cluster overdensity
  // this is the overdensity within which mass defined (i.e. \Delta)
  // set to -1.0 for virial radius, or 200 for M200 (rhocrit)
  float overden_id = -1.0; // 200 for delta=200 rho-c , -1 for delta=vir x rho-c
  int relation = 3; // concentration relation
  float rcutoff = 2.0;
	
  float Redshift = z;
  cosmic_t = cosm_model.cosmic_time(Redshift);
  cosmic_t0 = cosm_model.cosmic_time(0.0);
  E = cosm_model.Efact(Redshift);
	
  cluster nfwclus(Mvir, Redshift, overden_id, relation, cosm_model);
	
  //nfwclus.concentration(conc_norm, conc_mass_norm); // set halo concentration using M-c relation of Duffy et al (08)
  //M500 = nfwclus.get_mass_overden(500.0);// Msun
  //R500 = nfwclus.get_rad_overden(500.0);// (physical) Mpc
  //Rvir = nfwclus.get_radius();
	
  float cvir = conc_norm * c_vir_DK15_fast(z, Mvir*h);
  nfwclus.set_conc(cvir);
  M500 = nfwclus.get_mass_overden(500.0);// Msun
  R500 = nfwclus.get_rad_overden(500.0);// (physical) Mpc
  Rvir = nfwclus.get_radius();
	
  //cout << M500 << " " << R500 << " " << Rvir << endl;
	
  gas_model icm_mod(delta_rel, ad_index, eps_fb, eps_dm, fs_0, fs_alpha, pturbrad, delta_rel_zslope, delta_rel_n);
	
  icm_mod.calc_fs(M500, Omega_b/Omega_M, cosmic_t0, cosmic_t);
  icm_mod.evolve_pturb_norm(Redshift, rcutoff);
  icm_mod.set_nfw_params(Mvir, Rvir, nfwclus.get_conc(), nfwclus.get_rhoi(), R500);
  icm_mod.set_mgas_init(Omega_b/Omega_M);
  icm_mod.findxs();
	
  icm_mod.solve_gas_model(verbose, 1e-5);

  float npoly_mod, gamma_mod, nnt_mod;
  gamma_mod = gamma_mod0 * pow((1.0+Redshift),gamma_mod_zslope);
  npoly_mod = 1.0/(gamma_mod - 1.0 );
  nnt_mod = F.n_nt_mod;

  double r, pre;
  vector<double> pressure;

  for(int xi=0;xi<x.size();xi++){
    r = (double) x[xi]*R500;
    pre = icm_mod.returnP_mod2(r, R500, x_break, npoly_mod, nnt_mod, x_smooth) / 1e6; // keV/cm^3
    if(pre < 0.0) pre = 0.0;
    pressure.push_back(pre);
  }
	
  return pressure;
		
}

void free_FFTdata(){
		
  delete [] r;
  delete [] a;
  delete [] k;
  delete [] wsave;
	
}

void FFT_density_profile(double *output, double *bin, int nbin){

  double mu;
	
  double dlnr;
  double logrc=(logrmin+logrmax)/2.;
  double dlogr=(logrmax-logrmin)/double(N);
  dlnr=dlogr*log(10.);
  double nc=double(N+1)/2.;
	
  double logkc=log10(kr)-logrc;

  for(int i=0;i<N;i++){
    r[i]= pow(10.,logrc+(i-nc+1)*dlogr);
    k[i]= pow(10.,logkc+(i-nc+1)*dlogr);
  }

  int index[nbin+1];
  for(int ibin=0;ibin<nbin;ibin++){
    index[ibin] = (int)((log10(bin[ibin])-log10(k[0]))/dlogr);
  }
  double bin_edge = log10(bin[nbin-1])+(log10(bin[nbin-1])-log10(bin[nbin-2]));
  bin_edge = pow(10., bin_edge);
  index[nbin] = (int)((log10(bin_edge)-log10(k[0]))/dlogr);

  for(int ibin=0;ibin<nbin;ibin++){
    //fprintf(stdout,"%d: ", ibin);
    //fprintf(stdout,".");fflush(stdout);

    mu=+0.5;      
    q =-0.5;
      	
    for(int i=0;i<N;i++){
      double prof;
      if(r[i] > pow(10., xp[0]) && r[i] < pow(10., xp[Nspline-1])){
	double logx = log10(r[i]);
	splint(xp-1, yp-1, yp2-1, Nspline, logx, &prof);
	prof = pow(10., prof);
      }else{
	prof = 0;
      }
      a[i] = r[i] * r[i] * sqrt(M_PI/2) * prof;
    }

    //CALL FORTRAN SUBROUTINE, compute Fast Hankel Transform
    int ok;
    int clogical=0;
    ok=C2FLOGICAL(clogical);
	
    fhti_(& N, & mu, & q, & dlnr, & kr, & kropt, wsave, & ok);
	
    clogical=F2CLOGICAL(ok);
	
    fhtq_(& N, a, & dir, wsave);
    //END FORTRAN CALLING

    // cic
    int i1 = index[ibin];
    int i2 = index[ibin]+1;

    double f1 = a[i1];
    double f2 = a[i2];

    double b1 = k[i1];
    double b2 = k[i2];

    double res = (f2-f1)/(b2-b1)*(bin[ibin]-b1)+f1;
    output[ibin] = res/bin[ibin];           
  }
  //fprintf(stdout,"\n");

}
