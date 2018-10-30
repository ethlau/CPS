#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "cosmo.h"
#include "cluster.h"
#include "gas_model.h"
#include "xray.h"

#include "allvars.h"
#include "nrutil.h"
#include "proto.h"

#include "cfortran.h"

#define MAXBINS 500

double tarray[ntmax]; //keV
double zarray[nzmax]; //Solar unit
double rarray[nrmax]; 
double lambda_table[ntmax][nrmax];
double tres, zres, eres;

const double megapc = 3.0857e24; // in cm

using namespace std;

#define Nx 100
static double xp[Nx], yp[Nx], yp2[Nx];
static int Nspline;

std::vector<double> calc_Shaw_xray_emissivity_profile(cosmo cosm_model, float z, float Mvir, std::vector<float> x);

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

std::vector<double> calc_Flender_xray_emissivity_profile(cosmo cosm_model, float z, float Mvir, std::vector<float> x);

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
  double n_nt_mod; // fiducial : 0.80

  double clump0;
  double clump_zslope;
  double x_clump;
  double alpha_clump1;
  double alpha_clump2;

};

static struct Flender_param F;

double Mvir_to_mass ( cosmo cosm_model, float z, float Mvir, float overdensity );

double calc_Flender_xray_lum (cosmo cosm_model, float z, float Mvir, std::vector<float> x);

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

int main(int argc, char *argv[]){
	
  if ( argc!=3 ){
    fprintf(stderr,"usage: %s input OPT\n", argv[0]);
    fprintf(stderr,"Do you have apec_table.dat? YES:OPT=1, NO:OPT=0\n");
    exit(1);
  } 
	
  int opt_readin = atoi(argv[2]);
	
  FILE *fin;
  fin = fopen(argv[1], "r");
  if(fin == NULL){
    fprintf(stderr, "Cannot find %s\n", argv[1]);
    exit(1);
  }
	
  /* Cosmological parameters */
  float H0, Omega_M, Omega_b, wt, Omega_k, ns;
  fscanf(fin, "%f %f %f %f %f %f\n", &H0, &Omega_M, &Omega_b, &wt, &Omega_k, &ns);
	
  cout << H0 << " " << Omega_M << " " << Omega_b << " " << wt << " " << Omega_k << " " << ns << endl;
	
  /* cosmology class */
  cosmo cosm_model(H0, Omega_M, Omega_b, Omega_k, wt);
	
  char inputPk[256];
  fscanf(fin, "%s", inputPk);
	
  set_halo_conc_DK15(inputPk, Omega_M, Omega_b, wt, H0/100., ns);
	
  /* Flender et al paramaters */
  fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", 
	 &F.alpha0, &F.n_nt, &F.beta, &F.eps_f, &F.eps_DM, &F.f_star, &F.S_star, &F.A_C, 
	 &F.gamma_mod0, &F.gamma_mod_zslope, &F.x_break, &F.x_smooth, &F.n_nt_mod,

	 &F.clump0, &F.clump_zslope, &F.x_clump, &F.alpha_clump1, &F.alpha_clump2);

	
  int nzbin;
  float zmin, zmax;
  fscanf(fin, "%d %f %f\n", &nzbin, &zmin, &zmax);
  cout << nzbin << " " << zmin << " " << zmax << endl;
	
  int nmbin;
  float logMvir_min, logMvir_max;
  fscanf(fin, "%d %f %f\n", &nmbin, &logMvir_min, &logMvir_max);
  cout << nmbin << " " << logMvir_min << " " << logMvir_max << endl;
	
  float z, Mvir, x;
  std::vector<float> zlist, Mlist, xlist;
  /* set up for redshift, virial mass, radius bin */
  //redshift
  for(int i=0;i<nzbin;i++){
    z = (log(zmax)-log(zmin))/(nzbin-1)*(float)(1.*i) + log(zmin);
    z = exp(z);
    zlist.push_back(z);
  }
	
  //virial mass (Bryan & Norman) in [Msun], not in [Msun/h]
  for(int i=0;i<nmbin;i++){
    Mvir = (logMvir_max-logMvir_min)/(nmbin-1)*(float)(1.*i) + logMvir_min;
    Mvir = pow(10.0, Mvir);
    Mlist.push_back(Mvir);
  }
	
  //radius in unit of R500, x = r/R500
  float xmin = 1e-4;
  float xmax = 100.0;
  for(int i=0;i<Nx;i++){
    x = (log10(xmax)-log10(xmin))/Nx*(float)(1.*i+0.5) + log10(xmin);
    x = pow(10.0, x);
    xlist.push_back(x);
  }
	
  //set the output file name
  char outname[256];
  fscanf(fin, "%s\n", outname);
  cout << outname << endl;
  fclose(fin);
	
  set_lambda_table(tarray,rarray,lambda_table, opt_readin); 

  set_FFTlog_param();
	
  FILE *fp;	
  fp = fopen(outname, "wb");
  if(fp == NULL){
    printf("Cannot create %s.\n", outname);
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

  /*
  double total_flux = 0.0;
  for(int i=0;i<nzbin;i++){
    for(int j=0;j<nmbin;j++){
      double flux;
      flux = calc_Flender_xray_flux (cosm_model, zlist[i], Mlist[j], xlist);
      total_flux += flux;
    }
  }
  */

  /*
  FILE *fp_scaling;
  fp_scaling = fopen("lx_m.txt", "w");
  if(fp_scaling == NULL){
    printf("Cannot create lx_m.txt\n" );
    exit(1);
  }
  */	
	
  //fprintf(fp_scaling, "# M500 [Msun] z Lx [ergs/s]\n");
  for(int i=0;i<nzbin;i++){
    fprintf(stdout, ".");fflush(stdout);
    for(int j=0;j<nmbin;j++){
      //cout << zlist[i] << " " << Mlist[j] << endl;
      std::vector<double> emission;
      emission = calc_Flender_xray_emissivity_profile(cosm_model, zlist[i], Mlist[j], xlist); // ergs/s/cm^3/str
      double lum, M500;     
      lum = calc_Flender_xray_lum(cosm_model, zlist[i], Mlist[j], xlist); // ergs/s
      M500 = Mvir_to_mass(cosm_model, zlist[i], Mlist[j], 500.0);
      //fprintf(fp_scaling, "%e %e %e\n", M500, zlist[i], lum);
      double yp1 = 1.e31;
      double ypn = 1.e31;
      Nspline=0;
      for(int k=0;k<Nx;k++){
        if(emission[k] > 0) Nspline += 1;
      }
      if(Nspline > 2){
        int id =0;
        for(int k=0;k<Nx;k++){
          if(emission[k] > 0){
            xp[id] = log10(xlist[k]);
            yp[id] = log10(emission[k]);
            id += 1;
          }
        }
        double xout = pow(10., xp[Nspline-1]);
        spline(xp-1, yp-1, Nspline, yp1, ypn, yp2-1);
        FFT_density_profile(tab_yl_int, bin, Nx);	
        /*
        for(int k=0;k<Nx;k++){
            otab_yl_int[k] = x_l_integral(pow(10., tab_l_ls[k]), xout);
            if(isnan(tab_yl_int[k])){
                fprintf(stderr, "found NaN!\n");
                cout << "#z=" << zlist[i] << ",Mvir=" << Mlist[j] << " [Msun]" << endl;
                exit(1);
            }
        }
        */
        }else{
          for(int k=0;k<Nx;k++) tab_yl_int[k] = 0.0;
        }
        fwrite(tab_yl_int, sizeof(double), NPOINTS, fp);
      }
  }
  fclose(fp);
  //fclose(fp_scaling);
  fprintf(stdout, "DONE!\n");
	
  free_halo_conc_DK15();

  free_FFTdata();
	
  return 0;
}

double Mvir_to_mass ( cosmo cosm_model, float z, float Mvir, float overdensity ) {

  float mass;
  float Redshift = z;
  float overden_id = -1.0; // 200 for delta=200 rho-c , -1 for delta=vir x rho-c
  int relation = 3; // concentration relation
  float conc_norm = F.A_C;
  float h =cosm_model.get_H0()/100.0;

  cluster nfwclus(Mvir, Redshift, overden_id, relation, cosm_model);
  float cvir = conc_norm * c_vir_DK15_fast(z, Mvir*h);
  nfwclus.set_conc(cvir);
  mass = nfwclus.get_mass_overden(overdensity);// Msun
  return mass;
 
}

double calc_Flender_xray_lum (cosmo cosm_model, float z, float Mvir, std::vector<float> x){
	
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

  float clump0 = F.clump0;
  float clump_zslope = F.clump_zslope;
  float x_clump = F.x_clump;
  float alpha_clump1 = F.alpha_clump1;
  float alpha_clump2 = F.alpha_clump2;

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
	
  gas_model icm_mod(delta_rel, ad_index, eps_fb, eps_dm, fs_0, fs_alpha, pturbrad, delta_rel_zslope, delta_rel_n);
	
  icm_mod.calc_fs(M500, Omega_b/Omega_M, cosmic_t0, cosmic_t);
  icm_mod.evolve_pturb_norm(Redshift, rcutoff);
  icm_mod.set_nfw_params(Mvir, Rvir, nfwclus.get_conc(), nfwclus.get_rhoi(), R500);
  icm_mod.set_mgas_init(Omega_b/Omega_M);
  icm_mod.findxs();
	
  icm_mod.solve_gas_model(verbose, 1e-5);
	
  double Rmax = icm_mod.thermal_pressure_outer_rad()*R500;
  //double Yanl = icm_mod.calc_Y(R500, Rvir, Rmax);
	
  double r, emi, dlum;
  double luminosity = 0.0;
  double flux = 0.0;

  // distances in Mpc
  //double D_A = cosm_model.ang_diam(Redshift);
  double D_L = cosm_model.lum_dist(Redshift);

  float npoly_mod, gamma_mod;
  gamma_mod = gamma_mod0 * pow((1.0+Redshift),gamma_mod_zslope);
  npoly_mod = 1.0/(gamma_mod - 1.0 );

  double dvol[x.size()]; // radial shell vol in cm^3
  for(int xi=0;xi<x.size();xi++){
    r = (double) x[xi]*R500;
    if ( xi == 0 ) {
        dvol[xi] = 4.0*M_PI*pow(r*megapc, 3.0);
    } else {
        dvol[xi] = 4.0*M_PI*pow(r*megapc, 3.0) - dvol[xi-1];
    }
  }
  for(int xi=0;xi<x.size();xi++){
    // r in Mpc;
    r = (double) x[xi]*R500;
    if(r >= Rmax){
        emi = 0.0;
    } else{
        double ngas,pressure, kT, clump, clump1;
        pressure = icm_mod.returnPth_mod2(r, R500, x_break, npoly_mod, x_smooth); //keV cm^-3
        ngas = icm_mod.return_ngas_mod(r, R500, x_break, npoly_mod); //cm^-3
        kT = pressure/ngas; // keV
    
        clump1 = icm_mod.return_clumpf(r, R500, clump0, x_clump, alpha_clump1, alpha_clump2) - 1.0;
        clump1 *= pow(1.+Redshift, clump_zslope);
        clump = 1.0 + clump1;
        if (clump < 1.0) clump = 1.0;
        ngas *= sqrt(clump);

        emi = icm_mod.return_xray_emissivity(ngas, kT, Redshift); // ergs/s/cm^3
        //emi = icm_mod.calc_xray_emissivity(r, R500, Redshift); // ergs/s/cm^3
        dlum = emi * dvol[xi]; // ergs/s
    } 
    
    luminosity += dlum;		
  }
  flux = luminosity/(4.0*M_PI*D_L*D_L);	//ergs/s/cm^2
  return luminosity;	
}

std::vector<double> calc_Flender_xray_emissivity_profile(cosmo cosm_model, float z, float Mvir, std::vector<float> x){
	
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

  float clump0 = F.clump0;
  float clump_zslope = F.clump_zslope;
  float x_clump = F.x_clump;
  float alpha_clump1 = F.alpha_clump1;
  float alpha_clump2 = F.alpha_clump2;

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
  //printf("concentration_vir =%f\n",cvir);
  nfwclus.set_conc(cvir);
  M500 = nfwclus.get_mass_overden(500.0);// Msun
  R500 = nfwclus.get_rad_overden(500.0);// (physical) Mpc
  Rvir = nfwclus.get_radius();
	
  gas_model icm_mod(delta_rel, ad_index, eps_fb, eps_dm, fs_0, fs_alpha, pturbrad, delta_rel_zslope, delta_rel_n);
	
  icm_mod.calc_fs(M500, Omega_b/Omega_M, cosmic_t0, cosmic_t);
  icm_mod.evolve_pturb_norm(Redshift, rcutoff);
  icm_mod.set_nfw_params(Mvir, Rvir, nfwclus.get_conc(), nfwclus.get_rhoi(), R500);
  icm_mod.set_mgas_init(Omega_b/Omega_M);
  icm_mod.findxs();
	
  icm_mod.solve_gas_model(verbose, 1e-5);
	
  double Rmax = icm_mod.thermal_pressure_outer_rad()*R500;
  //double Yanl = icm_mod.calc_Y(R500, Rvir, Rmax);
	
  double r, emi;
  std::vector<double> emission;

  double fac = 4.0*M_PI; // in steradians

  float npoly_mod, gamma_mod;
  gamma_mod = gamma_mod0 * pow((1.0+Redshift),gamma_mod_zslope);
  npoly_mod = 1.0/(gamma_mod - 1.0 );

  for(int xi=0;xi<x.size();xi++){
    r = (double) x[xi]*R500; // in Mpc
    if(r >= Rmax){emi = 0.0;}
    else{
        double ngas,pressure, kT, clump, clump1;
        pressure = icm_mod.returnPth_mod2(r, R500, x_break, npoly_mod, x_smooth); //keV cm^-3
        ngas = icm_mod.return_ngas_mod(r, R500, x_break, npoly_mod); //cm^-3
        kT = pressure/ngas; // keV
    
        clump1 = icm_mod.return_clumpf(r, R500, clump0, x_clump, alpha_clump1, alpha_clump2) - 1.0;
        clump1 *= pow(1.+Redshift, clump_zslope);
        clump = 1.0 + clump1;
        if (clump < 1.0) clump = 1.0;
        ngas *= sqrt(clump);

        emi = icm_mod.return_xray_emissivity(ngas, kT, Redshift); // ergs/s/cm^3
        //printf("r, pressure, kT, ngas, emi = %f, %e, %e, %e, %e\n", r, pressure, kT, ngas, emi);
    } 		
    //if(emi < 1e-50){emi = 0.0;}
		
    emi = emi/fac/pow((1.0+Redshift),4.0); // ergs/s/cm^3/str
		
    //cout << r << " " << emi <<endl;`
		
    emission.push_back(emi);
  }
	
  return emission;
	
}


std::vector<double> calc_Shaw_xray_emissivity_profile(cosmo cosm_model, float z, float Mvir, std::vector<float> x){
	
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
	
  double Rmax = icm_mod.thermal_pressure_outer_rad()*R500;
  //double Yanl = icm_mod.calc_Y(R500, Rvir, Rmax);
	
  double r, emi;
  std::vector<double> emission;

  for(int xi=0;xi<x.size();xi++){
    r = (double) x[xi]*R500;
    if(r >= Rmax){emi = 0.0;}
    else{
        double pressure, ngas, kT;
        pressure = icm_mod.returnPth(r, R500); //keV cm^-3
        ngas = icm_mod.return_ngas(r); //cm^-3
        kT = pressure/ngas; // keV
        emi = icm_mod.return_xray_emissivity(ngas, kT, Redshift); // ergs/s/cm^3
    }
		
    //if(emi < 1e-50){emi = 0.0;}
		
    emi = emi/pow(1.+Redshift, 4.)/(4.0*M_PI); // ergs/s/cm^3/str
		
    //cout << r << " " << emi <<endl;
		
    emission.push_back(emi);
  }
	
  return emission;
	
}

/*
double x_l_integral_int_x(double x, double l_ls){
	double res;
	double x2px;
	if(x > pow(10., xp[0]) && x < pow(10., xp[Nspline-1])){
		double logx = log10(x);
		splint(xp-1, yp-1, yp2-1, Nspline, logx, &x2px);
		x2px = pow(10., x2px) *x*x;
	}else{
		return 0;
	}
	
	double kernel;
	double lx = l_ls*x;
	if(lx < 1e-2){
		kernel = 1.0-lx*lx/6.;
	}else{
		kernel = sin(lx)/(lx);
	}
	res =x2px*kernel;
	
	return res;
}
double x_l_integral(double l_ls, double xout){
	double res,abserr;
	size_t neval;
	
	int Nint = 10;
	//Romberg
	int i,j;
	double h[Nint];
	double s[Nint][Nint];
	
	res = 0;
	int NLoop = (int)((xout)/(2*M_PI/l_ls));
	
	if(NLoop < 2){
		double x_lo = 0.0;
		double x_hi = xout;
		
		for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}
		
		s[0][0] = 0.5*h[0]*(x_l_integral_int_x(x_hi, l_ls)+x_l_integral_int_x(x_lo, l_ls));
		for(i=2;i<=Nint;i++){
			s[i-1][0] = s[i-2][0];
			for(j=1;j<=pow(2.,i-2);j++){
				s[i-1][0] += h[i-2]*x_l_integral_int_x(x_lo+(2*j-1)*h[i-1], l_ls);
			}
			s[i-1][0] = 0.5*s[i-1][0];
		}
		
		for(i=2;i<=Nint;i++){
			for(j=2;j<=i;j++){
				s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
			}
		}
		
		res += s[Nint-1][Nint-1];
	}else{
		for(int iLoop=0;iLoop<NLoop;iLoop++){
			
			double x_lo = 0.0+(double)(iLoop+0)*(xout)/NLoop;
			double x_hi = 0.0+(double)(iLoop+1)*(xout)/NLoop;
			
			if(iLoop == NLoop-1){
				x_hi = xout;
			}
			
			for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}
			
			s[0][0] = 0.5*h[0]*(x_l_integral_int_x(x_hi, l_ls)+x_l_integral_int_x(x_lo, l_ls));
			for(i=2;i<=Nint;i++){
				s[i-1][0] = s[i-2][0];
				for(j=1;j<=pow(2.,i-2);j++){
					s[i-1][0] += h[i-2]*x_l_integral_int_x(x_lo+(2*j-1)*h[i-1], l_ls);
				}
				s[i-1][0] = 0.5*s[i-1][0];
			}
			
			for(i=2;i<=Nint;i++){
				for(j=2;j<=i;j++){
					s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
				}
			}
			
			res += s[Nint-1][Nint-1];
		}
	}
	return res;
}
*/

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
