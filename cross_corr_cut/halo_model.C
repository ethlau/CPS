#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "proto.h"
#include "allvars.h"
#include "nrutil.h"
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
#include <gsl/gsl_multiroots.h>

using namespace std;

static char TransTable[256];
static double xp[1000],yp[1000],yp2[1000];
static double tab_z[NPOINTS],tab_chi[NPOINTS],err_chi[NPOINTS];
static double scale_f[NPOINTS],GF[NPOINTS],GF2[NPOINTS];

static double tab_m[NPOINTS];
static double **tab_mf,**err_mf;
static double **tab_bi,**err_bi;

static double thre;
static double tab_R[NPOINTS],tab_dsdr[NPOINTS],tab_ds2dr2[NPOINTS],err_dsdr[NPOINTS],err_ds2dr2[NPOINTS];
static double tab_sig2[NPOINTS],err_sig2[NPOINTS];

static double **tab_r500,**err_r500;
static double **tab_r180m, **err_r180m;
static double **tab_rvir, **err_rvir;
static double **tab_r200, **err_r200;
static double **tab_mdeltatomvir, **err_mdeltatomvir;
static double **tab_cv,**err_cv;
static double **tab_nuM, **err_nuM;

static int WhichSpectrum,np,WhichWindow, OPT,WhichGF;
static int OPT_survey,OPT_fit;

static double bling, scale_Pk;

static double r_tophat,Rsmooth,Delta_c,fdelta_c,Const_MF;

static double AA,BB,CC;
static double B1,B2,B3,B4,B5;
static double nu, sigma, Omega_z;

static double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
static double G, Hubble;

static double Norm, InitTime, Redshift, a_exp;
static double Dplus; /* growth factor */

static double calc_z, z_source, multipole;


struct nfw_params{
	double wave_num;
	double con;
	double rad;
	double rho;
};

struct gas_params{
	double wave_num;
	double z;
};

struct fft_params{
	double wave_num;
	double con;
};

struct get_mvir{
	double delta;
	double rd;
	double redshift;
};

static double rho_mean, gal_norm;
//mainly based on arXiv:1106.3208 and arxiv:1412.5593

int mass_z_selection (double lmass, double z) {
    int value = 0;

    if (lmass < 12.0 + 2.5*(1.0+z)) value = 1;

    return value;
}


int main(int argc, char **argv) 
{
  int i,j,ip,iM,ir, idummy, nsnap;
  double dk,k,kmin,kmax,x;
  double TableMass[250],TabledN[250];
  double b_ling;
  char   outname[250];
	
  FILE *fred;
  FILE *fd;
  FILE *fdt;
	
  //input parameters======================================
  WhichSpectrum = 3;  // 1 Eisenstein-Hu, 2 Bond-Efstathiou, 3 CMBFAST table, 4 BBKS
  WhichGF = 2; // 1 apploximation formula, 2 caluculate linear density growth eq.	
  OPT = 6;     // 0 JK  1 ST  2 PS 3 Yahagi 4 Reed et al (astro-ph/0607150v4) 5 Bhattacharya et al (astro-ph/1005.2239) 6 Tinker et al 2008,2010
  OPT_fit = 2; // 1 smith et al 03, 2 takahashi et al 12
	
  if(argc!=10){
    fprintf(stderr,"usage:\n > halo_model outputname xl_file yl_file input_pk Om w0 h0 OPT_survey z_source\n");
    fprintf(stderr,"OPT_survey = 0 (delta_function), = 1 (for mimic HSC), =2 (for mimic LSST), = 3 (from COSMOS data), =4 (CFHTLenS)\n");
    return 1;
  }
	
  double hnu_kT = planck_h*(freq_obs*1e9)/(8.617342e-5)/T_CMB;
  double sz_fac = ((hnu_kT)*(exp(hnu_kT)+1)/(exp(hnu_kT)-1)-4.);
	
  cout << "g_nu = " << sz_fac << endl;
	
  sprintf(TransTable, "%s", argv[4]);
  Omega = atof(argv[5]);
  OmegaLambda = 1.0-Omega;
  w = atof(argv[6]);
  HubbleParam = atof(argv[7]);
	
  OPT_survey = atof(argv[8]);
  z_source   = atof(argv[9]);
	
  if(OPT_survey==3 || OPT_survey==4){gal_norm=gal_dist_norm(6.);}
  else{gal_norm=1;}
	
  fprintf(stdout, "Om = %4.3f, OL = %4.3f, w0 = %4.3f, h0 = %4.3f\n", Omega, OmegaLambda, w, HubbleParam);
  fprintf(stdout, "zs = %4.3f\n", z_source);
	
  if(WhichGF==2){growth_spline();}
	
  i=set_units();
  i=initialize_powerspectrum(WhichSpectrum);	
  stack_table_and_spline();
	
  stack_Pk_nonl_data_and_spline();
	
  sprintf(outname,"%s",argv[1]);
  fd=fopen(outname,"w");
  if(fd == NULL){
    fprintf(stderr, "can not make %s\n", outname);
    exit(1);
  }
	
  int Nbin=30;
  double lmax = 3.e4;
  double lmin = 10.0;
	
  double Mpc2cm = 3.0856*1e18*1e6;
  double yp1 = 1.e31;
  double ypn = 1.e31;
	
  fprintf(fd, "# ell cl_xk_1h cl_xk_2h cl_xx_1h cl_xx_2h cl_yk_1h cl_yk_2h cl_yy_1h cl_yy_2h cl_xy_1h cl_xy_2h cl_kk Nmode\n");
  for(i=0;i<Nbin;i++){
    double ell_bin_min = pow(10., log10(lmin) + (double)(i+0.0)*(log10(lmax)-log10(lmin))/Nbin);
    double ell_bin     = pow(10., log10(lmin) + (double)(i+0.5)*(log10(lmax)-log10(lmin))/Nbin);
    double ell_bin_max = pow(10., log10(lmin) + (double)(i+1.0)*(log10(lmax)-log10(lmin))/Nbin);
		
    FILE *fred1, *fred2;	
    fred1 = fopen(argv[2], "rb");
    if(fred1 == NULL){
      printf("you can not find %s.\n", argv[2]);
      exit(1);
    }
    fred2 = fopen(argv[3], "rb");
    if(fred2 == NULL){
      printf("you can not find %s.\n", argv[3]);
      exit(1);
    }
		
    int nzbin;
    float zmin, zmax;
    int nzbin_d;
    float zmin_d, zmax_d;
		
    fread(&nzbin, sizeof(int), 1, fred1);
    fread(&zmin, sizeof(float), 1, fred1);
    fread(&zmax, sizeof(float), 1, fred1);
		
    fread(&nzbin_d, sizeof(int), 1, fred2);
    fread(&zmin_d, sizeof(float), 1, fred2);
    fread(&zmax_d, sizeof(float), 1, fred2);
		
    if(nzbin_d != nzbin || fabs(zmin-zmin_d)>1e-5 || fabs(zmax-zmax_d)>1e-5){
      fprintf(stderr,"somthing wrong with input z tables!\n");
      cout << nzbin << " " << nzbin_d << endl;
      cout << zmin << " " << zmin_d << endl;
      cout << zmax << " " << zmax_d << endl;
      exit(1);
    }
		
    int nmbin;
    float logMvir_min;
    float logMvir_max;
    int nmbin_d;
    float logMvir_min_d;
    float logMvir_max_d;		
		
    fread(&nmbin, sizeof(int), 1, fred1);
    fread(&logMvir_min, sizeof(float), 1, fred1);
    fread(&logMvir_max, sizeof(float), 1, fred1);
		
    fread(&nmbin_d, sizeof(int), 1, fred2);
    fread(&logMvir_min_d, sizeof(float), 1, fred2);
    fread(&logMvir_max_d, sizeof(float), 1, fred2);
		
    if(nmbin_d != nmbin || fabs(logMvir_min-logMvir_min_d)>1e-5 || fabs(logMvir_max-logMvir_max_d)>1e-5){
      fprintf(stderr,"somthing wrong with input Mvir tables!\n");
      cout << nmbin << " " << nmbin_d << endl;
      cout << logMvir_min << " " << logMvir_min_d << endl;
      cout << logMvir_max << " " << logMvir_max_d << endl;
      exit(1);
    }
		
    int Nell;
    int Nell_d;
    fread(&Nell, sizeof(int), 1, fred1);
    fread(&Nell_d, sizeof(int), 1, fred2);
		
    if(Nell != Nell_d){
      fprintf(stderr,"somthing wrong with input ell tables!\n");
      cout << Nell << " " << Nell_d << endl;
      exit(1);
    }
		
    double tab_l_ls[Nell], tab_dummy[Nell];
    double tab_xl_int[Nell], tab_xl_int2[Nell];
    double tab_yl_int[Nell], tab_yl_int2[Nell];
		
    fread(tab_l_ls,  sizeof(double), Nell, fred1);
    fread(tab_dummy, sizeof(double), Nell, fred2);
		
    float dlnz = (log(zmax)-log(zmin))/(nzbin-1);
    float dlogm = (logMvir_max-logMvir_min)/(nmbin-1);
		
    double cl_xk_1=0.0, cl_xk_2=0.0;
    double cl_xx_1=0.0, cl_xx_2=0.0;
		
    double cl_yk_1=0.0, cl_yk_2=0.0;
    double cl_yy_1=0.0, cl_yy_2=0.0;
		
    double cl_xy_1=0.0, cl_xy_2=0.0;
		
    double zlist[nzbin];
		
    double cl_xk_1_int_z[nzbin], cl_xk_1_int_z2[nzbin];
    double cl_xk_2_int_z[nzbin], cl_xk_2_int_z2[nzbin];
    double cl_xx_1_int_z[nzbin], cl_xx_1_int_z2[nzbin];
    double cl_xx_2_int_z[nzbin], cl_xx_2_int_z2[nzbin];
		
    double cl_yk_1_int_z[nzbin], cl_yk_1_int_z2[nzbin];
    double cl_yk_2_int_z[nzbin], cl_yk_2_int_z2[nzbin];
    double cl_yy_1_int_z[nzbin], cl_yy_1_int_z2[nzbin];
    double cl_yy_2_int_z[nzbin], cl_yy_2_int_z2[nzbin];
		
    double cl_xy_1_int_z[nzbin], cl_xy_1_int_z2[nzbin];
    double cl_xy_2_int_z[nzbin], cl_xy_2_int_z2[nzbin];
		
    double mlist[nmbin];
		
    double cl_xk_1_int_m[nmbin], cl_xk_1_int_m2[nmbin];
    double cl_xx_1_int_m[nmbin], cl_xx_1_int_m2[nmbin];
    double cl_xx_2_int_m[nmbin], cl_xx_2_int_m2[nmbin];
		
    double cl_yk_1_int_m[nmbin], cl_yk_1_int_m2[nmbin];
    double cl_yy_1_int_m[nmbin], cl_yy_1_int_m2[nmbin];
    double cl_yy_2_int_m[nmbin], cl_yy_2_int_m2[nmbin];
		
    double cl_xy_1_int_m[nmbin], cl_xy_1_int_m2[nmbin];
    double cl_kk_2_int_m[nmbin], cl_kk_2_int_m2[nmbin];
		
    for(int iz=0;iz<nzbin;iz++){

      double zhere = exp(dlnz*(double)(1.*iz) + log(zmin));
      zlist[iz] = log(zhere);

      double covd = chi_fast(zhere);
      double dVdz = covd*covd * C*HubbleParam/H_z(zhere);
      double calc_k = ell_bin/covd;
      double gfac = (growth(1./(1+zhere))/growth(1.0));
      double Pk = gfac*gfac*PowerSpec(calc_k);
      double sigma_crit_inv = 4*PI*G0/C/C*covd*da_comb(zhere)/(1.+zhere); //(Msun/h)^-1

      for(int jm=0;jm<nmbin;jm++){

	double logMvir = dlogm*(double)(1.*jm) + logMvir_min;
	mlist[jm] = logMvir;

	double Mvir = pow(10., logMvir) * HubbleParam; // in the unit of Msun/h

	double m500=M500(zhere,Mvir); 
	double r500=R500_fast(zhere,Mvir);
	double ells = covd/(1+zhere)/r500; // = ell_500	
	double l_ls = ell_bin/ells;

	fread(tab_xl_int, sizeof(double), Nell, fred1);
	fread(tab_yl_int, sizeof(double), Nell, fred2);

	/*
	  for(int il=0;il<Nell;il++){
	  cout << il << " " << tab_xl_int[il] << " " << tab_yl_int[il] << endl;
	  }
	*/

	for(int il=0;il<Nell;il++){
	  tab_xl_int[il] = tab_xl_int[il] * Mpc2cm; // count/s/cm^2/arcsec^2/Mpc
	  //if(tab_xl_int[il] < 1e-70) tab_xl_int[il] = 1e-70;
					
	  tab_yl_int[il] = 1000*tab_yl_int[il]; // in the unit of eV/cc
	  tab_yl_int[il] = 4*M_PI* tab_yl_int[il] * sigma_T/m_elect * Mpc2cm; // Mpc^-1
	}
				
	double xl, yl;
	int index;
				
	index=-1;
	for (int il=0;il<Nell-1;il++) {
          if(tab_l_ls[il] < log10(l_ls) && log10(l_ls) < tab_l_ls[il+1]){
	    index = il;
          }
        }
				
	if(index < 0){
	  xl = 0;
	} else {
	  //xl = (log10(fabs(tab_xl_int[index+1]))-log10(fabs(tab_xl_int[index])))/(tab_l_ls[index+1]-tab_l_ls[index])
	  //  *(log10(l_ls)-tab_l_ls[index]) + log10(fabs(tab_xl_int[index]));
	  //xl = pow(10., xl) * 4*M_PI*(r500/HubbleParam)/ells/ells; // count/s/cm^2/arcsec^2
	  xl = (tab_xl_int[index+1]-tab_xl_int[index])/(tab_l_ls[index+1]-tab_l_ls[index])*(log10(l_ls)-tab_l_ls[index]) + tab_xl_int[index];
	  xl = xl * 4*M_PI*(r500/HubbleParam)/ells/ells; // count/s/cm^2/arcsec^2 
	}
	if(index < 0){
	  yl = 0;
	} else {
	  yl = (tab_yl_int[index+1]-tab_yl_int[index])/(tab_l_ls[index+1]-tab_l_ls[index])*(log10(l_ls)-tab_l_ls[index]) + tab_yl_int[index];
	  yl = yl * (r500/HubbleParam)/ells/ells; // dimensionless
	}
				
	double m200m = M_vir_to_M_delta(zhere, Mvir, 200.0);
	double mf = dndlogm_fast(log10(m200m), zhere);
	double b = halo_bias_fast(log10(m200m), zhere);
				
	double cnfw = c_nfw(zhere,Mvir);
	double rs   = r_s(zhere,Mvir);
	double rhos = rho_s(zhere,Mvir);
				
	struct nfw_params p = {calc_k,cnfw,(1.+zhere)*rs,rhos};
				
	double norm_m = 4*PI*rhos*rs*rs*rs*(log(1.+cnfw)-cnfw/(1.+cnfw));
	double Kappa = norm_m*U_nfw(Mvir, &p)*(1.+zhere)*(1.+zhere)/covd/covd;
	
        if ( logMvir <= 14.5 + 0.4*log10(zhere) ) {		
	    cl_xx_1_int_m[jm] = mf * xl * xl;
	    cl_xx_2_int_m[jm] = mf * b * xl;
	    cl_xk_1_int_m[jm] = mf * xl * Kappa;
				
	    cl_yy_1_int_m[jm] = mf * yl * yl;
	    cl_yy_2_int_m[jm] = mf * b * yl;
	    cl_yk_1_int_m[jm] = mf * yl * Kappa;
				
	    cl_kk_2_int_m[jm] = mf * b * Kappa;
				
	    cl_xy_1_int_m[jm] = mf * xl * yl;

        } else {
	    cl_xx_1_int_m[jm] = 0.0;
	    cl_xx_2_int_m[jm] = 0.0;
	    cl_xk_1_int_m[jm] = 0.0;
				
	    cl_yy_1_int_m[jm] = 0.0;
	    cl_yy_2_int_m[jm] = 0.0;
	    cl_yk_1_int_m[jm] = 0.0;
				
	    cl_kk_2_int_m[jm] = 0.0;
				
	    cl_xy_1_int_m[jm] = 0.0;
        }
      }
			
      spline(mlist-1, cl_yy_1_int_m-1, nmbin, yp1, ypn, cl_yy_1_int_m2-1);
      spline(mlist-1, cl_yy_2_int_m-1, nmbin, yp1, ypn, cl_yy_2_int_m2-1);
      spline(mlist-1, cl_yk_1_int_m-1, nmbin, yp1, ypn, cl_yk_1_int_m2-1);
			
      spline(mlist-1, cl_xx_1_int_m-1, nmbin, yp1, ypn, cl_xx_1_int_m2-1);
      spline(mlist-1, cl_xx_2_int_m-1, nmbin, yp1, ypn, cl_xx_2_int_m2-1);
      spline(mlist-1, cl_xk_1_int_m-1, nmbin, yp1, ypn, cl_xk_1_int_m2-1);
			
      spline(mlist-1, cl_kk_2_int_m-1, nmbin, yp1, ypn, cl_kk_2_int_m2-1);
			
      spline(mlist-1, cl_xy_1_int_m-1, nmbin, yp1, ypn, cl_xy_1_int_m2-1);
			
      double oneh_xx, twoh_xx, oneh_xk;
      double oneh_yy, twoh_yy, oneh_yk;
      double oneh_xy, twoh_kk;
			
      //double tab_spline_and_integral(int Nbin, double *xlist, double *ylist, double *zlist)
			
      oneh_xx = tab_spline_and_integral(nmbin, mlist, cl_xx_1_int_m, cl_xx_1_int_m2);
      twoh_xx = tab_spline_and_integral(nmbin, mlist, cl_xx_2_int_m, cl_xx_2_int_m2);
      oneh_xk = tab_spline_and_integral(nmbin, mlist, cl_xk_1_int_m, cl_xk_1_int_m2);
			
      oneh_yy = tab_spline_and_integral(nmbin, mlist, cl_yy_1_int_m, cl_yy_1_int_m2);
      twoh_yy = tab_spline_and_integral(nmbin, mlist, cl_yy_2_int_m, cl_yy_2_int_m2);
      oneh_yk = tab_spline_and_integral(nmbin, mlist, cl_yk_1_int_m, cl_yk_1_int_m2);			
			
      twoh_kk = tab_spline_and_integral(nmbin, mlist, cl_kk_2_int_m, cl_kk_2_int_m2);
      oneh_xy = tab_spline_and_integral(nmbin, mlist, cl_xy_1_int_m, cl_xy_1_int_m2);
			
      cl_xx_1_int_z[iz] = zhere * dVdz * oneh_xx;
      cl_xk_1_int_z[iz] = zhere * dVdz * oneh_xk * sigma_crit_inv;			
      cl_xx_2_int_z[iz] = zhere * dVdz * twoh_xx * twoh_xx * Pk;
      cl_xk_2_int_z[iz] = zhere * dVdz * twoh_xx * twoh_kk * Pk * sigma_crit_inv;
			
      cl_yy_1_int_z[iz] = zhere * dVdz * oneh_yy;
      cl_yk_1_int_z[iz] = zhere * dVdz * oneh_yk * sigma_crit_inv;			
      cl_yy_2_int_z[iz] = zhere * dVdz * twoh_yy * twoh_yy * Pk;
      cl_yk_2_int_z[iz] = zhere * dVdz * twoh_yy * twoh_kk * Pk * sigma_crit_inv;
			
      cl_xy_1_int_z[iz] = zhere * dVdz * oneh_xy;			
      cl_xy_2_int_z[iz] = zhere * dVdz * twoh_xx * twoh_yy * Pk;
			
    }
    fclose(fred1);
    fclose(fred2);
		
    spline(zlist-1, cl_xx_1_int_z-1, nzbin, yp1, ypn, cl_xx_1_int_z2-1);
    spline(zlist-1, cl_xx_2_int_z-1, nzbin, yp1, ypn, cl_xx_2_int_z2-1);
    spline(zlist-1, cl_xk_1_int_z-1, nzbin, yp1, ypn, cl_xk_1_int_z2-1);
    spline(zlist-1, cl_xk_2_int_z-1, nzbin, yp1, ypn, cl_xk_2_int_z2-1);
		
    spline(zlist-1, cl_yy_1_int_z-1, nzbin, yp1, ypn, cl_yy_1_int_z2-1);
    spline(zlist-1, cl_yy_2_int_z-1, nzbin, yp1, ypn, cl_yy_2_int_z2-1);
    spline(zlist-1, cl_yk_1_int_z-1, nzbin, yp1, ypn, cl_yk_1_int_z2-1);
    spline(zlist-1, cl_yk_2_int_z-1, nzbin, yp1, ypn, cl_yk_2_int_z2-1);
		
    spline(zlist-1, cl_xy_1_int_z-1, nzbin, yp1, ypn, cl_xy_1_int_z2-1);
    spline(zlist-1, cl_xy_2_int_z-1, nzbin, yp1, ypn, cl_xy_2_int_z2-1);
		
    cl_xx_1 = tab_spline_and_integral(nzbin, zlist, cl_xx_1_int_z, cl_xx_1_int_z2);
    cl_xx_2 = tab_spline_and_integral(nzbin, zlist, cl_xx_2_int_z, cl_xx_2_int_z2);
    cl_xk_1 = tab_spline_and_integral(nzbin, zlist, cl_xk_1_int_z, cl_xk_1_int_z2);
    cl_xk_2 = tab_spline_and_integral(nzbin, zlist, cl_xk_2_int_z, cl_xk_2_int_z2);
		
    cl_yy_1 = tab_spline_and_integral(nzbin, zlist, cl_yy_1_int_z, cl_yy_1_int_z2);
    cl_yy_2 = tab_spline_and_integral(nzbin, zlist, cl_yy_2_int_z, cl_yy_2_int_z2);
    cl_yk_1 = tab_spline_and_integral(nzbin, zlist, cl_yk_1_int_z, cl_yk_1_int_z2);
    cl_yk_2 = tab_spline_and_integral(nzbin, zlist, cl_yk_2_int_z, cl_yk_2_int_z2);
		
    cl_xy_1 = tab_spline_and_integral(nzbin, zlist, cl_xy_1_int_z, cl_xy_1_int_z2);
    cl_xy_2 = tab_spline_and_integral(nzbin, zlist, cl_xy_2_int_z, cl_xy_2_int_z2);
		
    double cl_kk = Pkappa(ell_bin, 2);
    double Nmode = 2*ell_bin*(ell_bin_max-ell_bin_min);
		
    fprintf(stdout,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n", ell_bin, 
	    cl_xk_1, cl_xk_2, 
	    cl_xx_1, cl_xx_2, 
	    cl_yk_1, cl_yk_2, 
	    cl_yy_1, cl_yy_2, 
	    cl_xy_1, cl_xy_2, 
	    cl_kk, 
	    Nmode);
    fprintf(fd,"%e %e %e %e %e %e %e %e %e %e %e %e %e\n", ell_bin, 
	    cl_xk_1, cl_xk_2, 
	    cl_xx_1, cl_xx_2, 
	    cl_yk_1, cl_yk_2, 
	    cl_yy_1, cl_yy_2, 
	    cl_xy_1, cl_xy_2, 
	    cl_kk, 
	    Nmode);
  }
  fclose(fd);
	
  free_dmatrix(tab_mf,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_mf,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_bi,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_bi,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_r500,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_r500,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_rvir,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_rvir,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_mdeltatomvir,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_mdeltatomvir,1,NPOINTS,1,NPOINTS);
	
  free_dmatrix(tab_nuM,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_nuM,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_cv,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_cv,1,NPOINTS,1,NPOINTS);
	
  return 0;
	
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
	Sigma8 = sqrt(TopHatSigma2(8.0));
	fprintf(stdout,"sigma8=%g \n",Sigma8);
	
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
	return Norm*pow(k,ns)/ pow(1+pow(AA*k+pow(BB*k,1.5)+CC*CC*k*k,nu),2/nu);
}




//==================================================
inline double PowerSpec_BBKS(double k)
{
	return Norm*pow(k,ns)* pow(log(1.0+B1*k)/(B1*k),2)/ 
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
	return Norm*transfunc_cmbfast(k);
}


//==================================================
double transfunc_cmbfast(double k)
{
	int i;
	double lk;
	double pow_index;
	
	//k *= 1000.0; /* convert to h/Mpc */
	
	if(k < 1e-10){
		return 0.0;
	}
	
	lk=log10(k);
	
	if(lk < xp[0]){
		double dummy = (yp[1]-yp[0])/(xp[1]-xp[0])*(lk-xp[0])+yp[0];
		return pow(10.,dummy);
	}
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
	char tfunc_table[80];
	int  errorFlag=0;
	FILE *fd;
	
	fprintf(stdout,"Reading %s .\n",TransTable);
	iline=0;
	if((fd=fopen(TransTable,"r"))!=NULL){
		while(!feof(fd)){
			fscanf(fd,"%lf %lf \n",&xp[iline],&yp[iline]);
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
	double kr,kr3,kr2,x;
	
	kr=r_tophat*k;
	kr2=kr*kr;
	kr3=kr2*kr;
	
	if(kr<1e-8) return 0;
	
	double window=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
	x=4*PI*k*k*window*window*PowerSpec(k)/pow(2*PI,3);
	
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
		return growth_for_any_w(a);
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
	
	res=dweight(x, rsphere)*PowerSpec(x);
	
	return res;
}

inline double var2(double x, double rsphere)
{
	double res,pspec,xk;
	
	//xk=x/1000.0; // kpc -> Mpc scaling 
	//res=weight(x, rsphere)*PowerSpec(xk) * (1000.0/Norm); //renormalize for mf code
	
	res= weight(x, rsphere)*PowerSpec(x);
	
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
	
	return res/pow(2*PI, 3);
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
			double Delta=delta_h;
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

double RungeKutta(double a_in, double a){
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
		GF[i] = RungeKutta(1.0/(1.0+1088.2), scale_f[i]);
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
	
	splint(scale_f-1, GF-1, GF2-1, NPOINTS, la, &pow_index);
	return pow(10.0,pow_index);
}

void stack_table_and_spline(){
	int i,j;
	double yp1,ypn;
	double dlogz = (log10(zhalo_max)-(-3.))/(NPOINTS-1);
	double dlogm = (log10(Mhalo_max)-7.5)/(NPOINTS-1);
	
	yp1 = 1.e31;
	ypn = 1.e31;
	
	
	for(i=0;i<NPOINTS;i++){
		tab_z[i] = i*dlogz + (-3.);
		tab_m[i] = i*dlogm + 7.5;
	}	
	
	for(i=0;i<NPOINTS;i++){
		tab_chi[i] = chi(pow(10,tab_z[i]));
		tab_chi[i] = log10(tab_chi[i]);
	}
	
	printf("stock chi(z) data...\n");
	spline(tab_z-1, tab_chi-1, NPOINTS, yp1, ypn, err_chi-1);
	printf("set spline chi(z)...\n");
	
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
	printf("set spline c200c(z, M)...\n");
	
	tab_mf=dmatrix(1,NPOINTS,1,NPOINTS);
	err_mf=dmatrix(1,NPOINTS,1,NPOINTS);
	tab_bi=dmatrix(1,NPOINTS,1,NPOINTS);
	err_bi=dmatrix(1,NPOINTS,1,NPOINTS);
	tab_r500=dmatrix(1,NPOINTS,1,NPOINTS);
	err_r500=dmatrix(1,NPOINTS,1,NPOINTS);
	tab_rvir=dmatrix(1,NPOINTS,1,NPOINTS);
	err_rvir=dmatrix(1,NPOINTS,1,NPOINTS);
	tab_r180m=dmatrix(1,NPOINTS,1,NPOINTS);
	err_r180m=dmatrix(1,NPOINTS,1,NPOINTS);
	tab_r200=dmatrix(1,NPOINTS,1,NPOINTS);
	err_r200=dmatrix(1,NPOINTS,1,NPOINTS);
	tab_mdeltatomvir=dmatrix(1,NPOINTS,1,NPOINTS);
	err_mdeltatomvir=dmatrix(1,NPOINTS,1,NPOINTS);
	
	for(i=0;i<NPOINTS;i++){
		InitTime=1.0/(1.0+pow(10,tab_z[i])); // expansion parameter at which MF is computed
		calc_z=1.0/InitTime -1.0;
		//cout << calc_z << endl;
		Delta_c=1.68647;
		a_exp=1.0/(1.0+calc_z);
		
		Omega_z = Omega*pow(1.0+calc_z,3.0)/(Omega*pow(1.0+calc_z,3.0) + OmegaLambda);
		fdelta_c = delta_c_func();
		//double b_ling=GrowthFactor(1.,InitTime)/(InitTime/1.0);
		//Rsmooth = 8.0; //here use Mpc/h
		//Const_MF = Sigma8*Sigma8/unnsigma(Rsmooth)*b_ling*b_ling;
		
		double b_ling=GrowthFactor(1.,InitTime);
		Const_MF = b_ling*b_ling;
		
		for(j=0;j<NPOINTS;j++){
			tab_mf[i+1][j+1] = dndlogm(tab_m[j]);
			tab_bi[i+1][j+1] = halo_bias(tab_m[j]);
			tab_rvir[i+1][j+1] = r_vir(pow(10,tab_z[i]), pow(10.,tab_m[j]));
			tab_r500[i+1][j+1] = R500(pow(10,tab_z[i]), pow(10.,tab_m[j]));
			tab_r200[i+1][j+1] = R200(pow(10,tab_z[i]), pow(10.,tab_m[j]));
			tab_r180m[i+1][j+1] = R180m(pow(10,tab_z[i]), pow(10.,tab_m[j]));
			
			tab_mdeltatomvir[i+1][j+1] = M_delta_to_M_vir(pow(10,tab_z[i]), pow(10.,tab_m[j]), delta_h);
			if(tab_mf[i+1][j+1] < 1e-100){
				printf("z=%e M=%e\n",pow(10,tab_z[i]),pow(10,tab_m[j]));
				printf("dndlogm=%e, bh=%e\n",tab_mf[i+1][j+1],tab_bi[i+1][j+1]);
				//printf("rvir=%e, r500=%e\n",tab_rvir[i+1][j+1],tab_r500[i+1][j+1]);
			}
		}
	}
	//exit(1);
	
	for(i=1;i<=NPOINTS;i++){
		for(j=1;j<=NPOINTS;j++){
			tab_mf[i][j] = log10(tab_mf[i][j]);
			tab_bi[i][j] = log10(tab_bi[i][j]);
			tab_r500[i][j] = log10(tab_r500[i][j]);
			tab_rvir[i][j] = log10(tab_rvir[i][j]);
			tab_r180m[i][j] = log10(tab_r180m[i][j]);
			tab_r200[i][j] = log10(tab_r200[i][j]);
			tab_mdeltatomvir[i][j] = log10(tab_mdeltatomvir[i][j]);
		}
	}
	
	printf("stock dn/dlogm data...\n");
	
	splie2(tab_z-1, tab_m-1,tab_mf,NPOINTS,NPOINTS,err_mf);
	printf("set spline dn/dlogm...\n");
	splie2(tab_z-1, tab_m-1,tab_bi,NPOINTS,NPOINTS,err_bi);
	printf("set spline halo bias...\n");
	splie2(tab_z-1, tab_m-1,tab_rvir,NPOINTS,NPOINTS,err_rvir);
	printf("set spline rvir...\n");
	splie2(tab_z-1, tab_m-1,tab_r500,NPOINTS,NPOINTS,err_r500);
	printf("set spline r500...\n");
	splie2(tab_z-1, tab_m-1, tab_r180m,NPOINTS,NPOINTS,err_r180m);
	printf("set spline r180m...\n");
	splie2(tab_z-1, tab_m-1, tab_r200,NPOINTS,NPOINTS,err_r200);
	printf("set spline r200...\n");
	splie2(tab_z-1, tab_m-1, tab_mdeltatomvir,NPOINTS,NPOINTS,err_mdeltatomvir);
	printf("set spline m200->mvir..\n");
	
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

double chi_fast(double z){
	double lz;
	double pow_index;
	
	lz=log10(z);
	
	if(lz < tab_z[0])
	return 1e-30; /* usually should not be executed */
	if(lz > tab_z[NPOINTS-1]){
		return chi(z);
	}
	
	splint(tab_z-1, tab_chi-1, err_chi-1, NPOINTS, lz, &pow_index);
	//printf("OK.spline chi(z)\n");
	return pow(10.0,pow_index);
}
double dndlogm_fast(double logm, double z){
	double lz = log10(z);
	double pow_index;
	
	splin2(tab_z-1, tab_m-1, tab_mf, err_mf, NPOINTS, NPOINTS, lz, logm, &pow_index);
	//printf("OK.spline dn/dlogm\n");
	if(isnan(pow_index)){
		printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
		exit(1);
	}
	return pow(10.0,pow_index);
}

double halo_bias_fast(double logm, double z){
	double lz = log10(z);
	double pow_index;
	
	splin2(tab_z-1, tab_m-1, tab_bi, err_bi, NPOINTS, NPOINTS, lz, logm, &pow_index);
	//printf("OK.spline dn/dlogm\n");
	if(isnan(pow_index))printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
	return pow(10.0,pow_index);
}

double R500_fast(double z, double M){
	double lz = log10(z);
	double logm = log10(M);
	
	double pow_index;
	
	if(lz > tab_z[NPOINTS-1] || lz < tab_z[0] || logm > tab_m[NPOINTS-1] || logm < tab_m[0]){
		return R500(z,M);
	}else{
		splin2(tab_z-1, tab_m-1, tab_r500, err_r500, NPOINTS, NPOINTS, lz, logm, &pow_index);
		//printf("OK.spline dn/dlogm\n");
		if(isnan(pow_index))printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
		return pow(10.0,pow_index);
	}
	
}

double R200_fast(double z, double M){
	double lz = log10(z);
	double logm = log10(M);
	
	double pow_index;
	
	if(lz > tab_z[NPOINTS-1] || lz < tab_z[0] || logm > tab_m[NPOINTS-1] || logm < tab_m[0]){
		return R200(z,M);
	}else{
		splin2(tab_z-1, tab_m-1, tab_r200, err_r200, NPOINTS, NPOINTS, lz, logm, &pow_index);
		//printf("OK.spline dn/dlogm\n");
		if(isnan(pow_index))printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
		return pow(10.0,pow_index);
	}
	
}

double R180m_fast(double z, double M){
	double lz = log10(z);
	double logm = log10(M);
	
	double pow_index;
	
	if(lz > tab_z[NPOINTS-1] || lz < tab_z[0] || logm > tab_m[NPOINTS-1] || logm < tab_m[0]){
		return R180m(z,M);
	}else{
		splin2(tab_z-1, tab_m-1, tab_r180m, err_r180m, NPOINTS, NPOINTS, lz, logm, &pow_index);
		//printf("OK.spline dn/dlogm\n");
		if(isnan(pow_index))printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
		return pow(10.0,pow_index);
	}
	
}

double rvir_fast(double z, double M){
	double lz = log10(z);
	double logm = log10(M);
	
	double pow_index;
	
	if(lz > tab_z[NPOINTS-1] || lz < tab_z[0] || logm > tab_m[NPOINTS-1] || logm < tab_m[0]){
		return r_vir(z,M);
	}else{
		splin2(tab_z-1, tab_m-1, tab_rvir, err_rvir, NPOINTS, NPOINTS, lz, logm, &pow_index);
		//printf("OK.spline dn/dlogm\n");
		if(isnan(pow_index))printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
		return pow(10.0,pow_index);
	}
	
}

double halo_bias(double logm){ //based on peak-background split
	double delta=1.68647;
	double sqrt_two_over_pi=0.79788456;
	double sqrt_two=1.414213562;
	double fdeltac;
	double rho_crit=2.7754e11;
	double rm,rmass,sig,ans,fract;
	double r, res;
	
	fdeltac=delta_c_func();
	
	rm = pow(10,logm);
	rmass = rm/rho_crit/Omega;
	
	sig = sigma_m(rmass, &r);
	
	double Nu = fdeltac/sig;
	
	//For press-schechter
	//res = 1.+ (nu*nu-1.)/fdeltac;
	
	if(OPT==6){
		//Tinker et al 2010 arXiv:1001.3162
		Nu = 1.686/sig;
		double a1,a2,b1,b2,c1,c2;
		double Delta=delta_h;
		double y = log10(Delta);
		a1 = 1.0+0.24*y*exp(-pow(4./y,4.));
		a2 = 0.44*y-0.88;
		b1 = 0.183;
		b2 = 1.5;
		c1 = 0.019+0.107*y+0.19*exp(-pow(4./y,4.));
		c2 = 2.4;
		res = 1.0-a1*(pow(Nu,a2))/(pow(Nu,a2)+pow(1.686,a2))+b1*pow(Nu,b2)+c1*pow(Nu,c2);
	}
	else if(OPT==5){
		
		Nu = fdeltac/sig;
		double A,B,p,q;
		A = 0.333/pow(1+calc_z,0.11);
		B = 0.788/pow(1+calc_z,0.01);
		p = 0.807;
		q = 1.795;
		
		double X=B*Nu*Nu;
		
		res = 1.0 + 2*X/fdeltac*(0.5 - 1/X*(q/2+(-p+q/2)*pow(X,-p))/(1+pow(X,-p)));
		
	}else{
		//For sheth-tormen
		Nu = fdeltac/sig;
		double a = 0.75,p=0.3;
		res = 1.0 + (a*Nu*Nu-1.)/fdeltac + (2.*p)/fdeltac/(1.+pow(a*Nu*Nu,p));
	}
	
	return res;
}

double H_z(double z){
	double res;
	
	if(w == -1.0){
		res=100*HubbleParam*sqrt(Omega*pow((1+z),3)+OmegaLambda);
	}else{
		res=100*HubbleParam*sqrt(Omega*pow((1+z),3)+OmegaLambda*pow(1+z,3*(1+w)));
	}
	
	return res;
}

double H_z_1(double z,void *params){
	double res;
	res=1/H_z(z);
	return res;
}

double integral_chi(double lnz, void *params){
	double z = exp(lnz);
	double res = z/H_z(z);
	
	return res;
}

double chi(double z){
	
	double result,abserr,params=0;
	size_t neval;
	gsl_function F;
	
	/*	
	F.function=&H_z_1;
	F.params=&params;
	if(z>1e-10){
	gsl_integration_qng(&F,0,z,0,1e-6,&result, &abserr,&neval);
	}else{
	result=0;
	}
	*/
	
	F.function=&integral_chi;
	F.params=&params;
	if(z>1e-10){
		gsl_integration_qng(&F,-10.,log(z),0,1e-6,&result, &abserr,&neval);
	}else{
		result=0;
	}
	
	return C*result*HubbleParam; //the unit of chi should be Mpc/h
}

double rho_nfw(double r,double z,double M){
	double x=r/r_s(r,M);
	return rho_s(z,M)/(x*pow((1+x),2));
}

double rho_s(double z, double M){
	//M = Mvir
	double cvir = c_nfw(z,M);
	double rs = r_s(z,M);
	
	return M/(log(1.+cvir)-cvir/(1.+cvir))/(4*PI*rs*rs*rs);
}

double r_s(double z,double M){
	return r_vir(z,M)/c_nfw(z,M);
}

double c_nfw(double z, double M){
	/*
	double logM,M1;
	double res;
	
	logM=log10(M)-14.0;
	M1=pow(10.0,logM);	
	res = 5.72/pow(1+z,0.71)*pow(M1,-0.081); // M (M_sun/h) cf.duffy et al (2008)	
	
	//see http://arxiv.org/pdf/1402.7073v2.pdf
	//
	//res = 0.537 + 0.488 * exp(-0.718*pow(z, 1.08)) + (-0.097+0.024*z)*log10(M/2e12);
	//res = pow(10., res);
	//if(res < 1.){res = 1.;}
	
	return res;
	*/
	return c_vir_DK15_fast(z, M);
}

double delta_c(double z,double M){
	return delta_v(z)/3.0*pow(c_nfw(z,M),3)/(log(1.0+c_nfw(z,M))-c_nfw(z,M)/(1.0+c_nfw(z,M)));
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

double M500(double z, double M){//See, e.g., Appendix D15 of Komatsu et al., arXiv:1001.4538
	double r500=R500_fast(z,M);
	double rs=r_s(z,M);
	double cnfw=c_nfw(z,M);
	
	double x=r500/rs;
	
	return M*(log(1+x)-x/(1+x))/(log(1+cnfw)-cnfw/(1+cnfw));
}

double R500(double z, double M){
	//printf("calc R500\n");
	
	//The definition of R500 (as follows in http://arxiv.org/pdf/0910.1234v1.pdf)
	// M500 = 4*PI/3 * 500 * rho_crit(z) * R500^3 = M_NFW (<R500) 
	
	int status;
	int iter=0, max_iter=500;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r=0;
	double rho_crit=2.7754e11; //M_sun/h/(Mpc/h)^3
	double cnfw = c_nfw(z,M);
	double rs = r_vir(z,M)/cnfw;
	double mc = log(1+cnfw)-cnfw/(1+cnfw);
	double param = 4*PI/3*500*rho_crit*(Omega*pow(1.+z,3)+OmegaLambda)*rs*rs*rs*mc/M;
	
	double x_lo=1e-5,x_hi=100*c_nfw(z,M);
	gsl_function F;
	F.function=&calc_R500_func;
	F.params=&param;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	
	gsl_root_fsolver_set(s,&F,x_lo,x_hi);
	
	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
		//if (status == GSL_SUCCESS) printf ("Converged:\n"); 
		//printf("%5d [%.7f, %.7f] %.7f %.7f\n",
		//		iter, x_lo, x_hi, r, x_hi-x_lo);
	} while (status == GSL_CONTINUE /*&& iter < max_iter*/);
	gsl_root_fsolver_free(s); 
	
	//printf("done R500\n");
	return r*rs;
}

double calc_R500_func(double x, void *params){
	double p = *(double*) params;
	double res = (log(1+x)-x/(1+x))-p*x*x*x;
	
	return res;
}

double M200(double z, double M){//See, e.g., Appendix D15 of Komatsu et al., arXiv:1001.4538
	double r500=R200_fast(z,M);
	double rs=r_s(z,M);
	double cnfw=c_nfw(z,M);
	
	double x=r500/rs;
	
	return M*(log(1+x)-x/(1+x))/(log(1+cnfw)-cnfw/(1+cnfw));
}

double R200(double z, double M){
	//printf("calc R500\n");
	
	//The definition of R200 (as follows in http://arxiv.org/pdf/1412.5593.pdf)
	// M200 = 4*PI/3 * 200 * rho_crit(z) * R200^3 = M_NFW (<R200) 
	
	int status;
	int iter=0, max_iter=500;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r=0;
	double rho_crit=2.7754e11; //M_sun/h/(Mpc/h)^3
	double cnfw = c_nfw(z,M);
	double rs = r_vir(z,M)/cnfw;
	double mc = log(1+cnfw)-cnfw/(1+cnfw);
	double param = 4*PI/3*200*rho_crit*(Omega*pow(1.+z,3)+OmegaLambda)*rs*rs*rs*mc/M;
	
	double x_lo=1e-5,x_hi=100*c_nfw(z,M);
	gsl_function F;
	F.function=&calc_R500_func;
	F.params=&param;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	
	gsl_root_fsolver_set(s,&F,x_lo,x_hi);
	
	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
		//if (status == GSL_SUCCESS) printf ("Converged:\n"); 
		//printf("%5d [%.7f, %.7f] %.7f %.7f\n",
		//		iter, x_lo, x_hi, r, x_hi-x_lo);
	} while (status == GSL_CONTINUE /*&& iter < max_iter*/);
	gsl_root_fsolver_free(s); 
	
	//printf("done R500\n");
	return r*rs;
}

double R180m(double z, double M){
	
	//The definition of R180m (as follows in 0205468v3.pdf)
	// M180m = 4*PI/3 * 180 * rho_m(z) * R180m^3 = M_NFW (<R180m) 
	
	int status;
	int iter=0, max_iter=500;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r=0;
	double rho_crit=2.7754e11; //M_sun/h/(Mpc/h)^3
	double cnfw = c_nfw(z,M);
	double rs = r_vir(z,M)/cnfw;
	double mc = log(1+cnfw)-cnfw/(1+cnfw);
	double param = 4*PI/3*180*rho_crit*(Omega*pow(1.+z,3))*rs*rs*rs*mc/M;
	
	double x_lo=1e-5,x_hi=100*c_nfw(z,M);
	gsl_function F;
	F.function=&calc_R500_func;
	F.params=&param;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	
	gsl_root_fsolver_set(s,&F,x_lo,x_hi);
	
	do {
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
		//if (status == GSL_SUCCESS) printf ("Converged:\n"); 
		//printf("%5d [%.7f, %.7f] %.7f %.7f\n",
		//		iter, x_lo, x_hi, r, x_hi-x_lo);
	} while (status == GSL_CONTINUE /*&& iter < max_iter*/);
	gsl_root_fsolver_free(s); 
	
	return r*rs;
}

double M180m(double z, double M){
	double r500=R180m_fast(z,M);
	double rs=r_s(z,M);
	double cnfw=c_nfw(z,M);
	
	double x=r500/rs;
	
	return M*(log(1+x)-x/(1+x))/(log(1+cnfw)-cnfw/(1+cnfw));
}

/*
double onehalo_int_gas(double logm, void *params){
double rm = pow(10,logm);

struct gas_params p = *(struct gas_params *)params;
double ell = p.wave_num;
double z = p.z;

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(z, rm);
logm = log10(m180m);
}
double mf = dndlogm_fast(logm, z);
double u  = y_l(ell, z, rm);

double res = mf * u * u;
return res;
}
double twohalo_int_gas(double logm, void *params){
double rm = pow(10,logm);

struct gas_params p = *(struct gas_params *)params;
double ell = p.wave_num;
double z = p.z;

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(z, rm);
logm = log10(m180m);
}

double mf = dndlogm_fast(logm, z);
double u  = y_l(ell, z, rm);
double b = halo_bias_fast(logm, z);

double res = mf*b*u;

return res;
}

double P_1halo_gas(double k, void *params){
double res,abserr;
size_t neval;
double x_lo = log10(5.0e11);
double x_hi = log10(5.0e15);

double z = *(double *)params;
struct gas_params p_int = {k, z};

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(onehalo_int_gas(x_hi,&p_int)+onehalo_int_gas(x_lo,&p_int));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*onehalo_int_gas(x_lo+(2*j-1)*h[i-1],&p_int);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}

double P_2halo_gas(double k, void *params){
double res,abserr;
double x_lo = log10(5.0e11);
double x_hi = log10(5.0e15);

double z = *(double *)params;
double dA = chi_fast(z);
double gfac=(growth(1./(1+z))/growth(1.0));
struct gas_params p_int = {k,z};

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(twohalo_int_gas(x_hi,&p_int)+twohalo_int_gas(x_lo,&p_int));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*twohalo_int_gas(x_lo+(2*j-1)*h[i-1],&p_int);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res*res*gfac*gfac*PowerSpec(k/dA);
}


double Cl_1halo_yy_int(double lnx, void *params){
double res;
double ell = *(double *)params;
double z = exp(lnx); //lnx = log(z)
double covd = chi_fast(z);

if(z > 1e-4){
res = z*covd*covd*C*HubbleParam/H_z(z)*P_1halo_gas_fast(ell, z); //dominate 1-halo term for ell > 300
}else{
res=1e-30;
}
//printf("%g %g %g\n",ell,z,res);

return res;

}

double Cl_2halo_yy_int(double lnx, void *params){
double res;
double ell = *(double *)params;
double z = exp(lnx); //lnx = log(z)
double covd = chi_fast(z);

if(z > 1e-4){
res = z*covd*covd*C*HubbleParam/H_z(z)*P_2halo_gas_fast(ell, z);
}else{
res=1e-30;
}
//printf("%g %g %g\n",ell,z,res);

return res;

}

double Cl_1halo_yy(double ell, double zmax){
double res,abserr;
size_t neval;
double x_lo = log(5e-3);
double x_hi = log(zmax);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(Cl_1halo_yy_int(x_hi, &ell)+Cl_1halo_yy_int(x_lo, &ell));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*Cl_1halo_yy_int(x_lo+(2*j-1)*h[i-1], &ell);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}

double Cl_2halo_yy(double ell, double zmax){
double res,abserr;
size_t neval;
double x_lo = log(5e-3);
double x_hi = log(zmax);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(Cl_2halo_yy_int(x_hi,&ell)+Cl_2halo_yy_int(x_lo,&ell));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*Cl_2halo_yy_int(x_lo+(2*j-1)*h[i-1],&ell);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}

double P_1halo_gas_fast(double ell, double z){
double lnz = log10(z);
double lnell = log10(ell); 

double pow_index;

if(lnz > tab_z[NPOINTS-1] || lnz < tab_z[0] || ell > pow(10., tab_ell[NPOINTS-1]) || ell < pow(10., tab_ell[0])){
return P_1halo_gas(ell, &z);
}else{
splin2(tab_ell-1, tab_z-1, tab_Cl_yy_1, err_Cl_yy_1, NPOINTS, NPOINTS, lnell, lnz, &pow_index);
return pow(10.0,pow_index);
}
}

double P_2halo_gas_fast(double ell, double z){
double lnz = log10(z);
double lnell = log10(ell); 

double pow_index;

if(lnz > tab_z[NPOINTS-1] || lnz < tab_z[0] || ell > pow(10., tab_ell[NPOINTS-1]) || ell < pow(10., tab_ell[0])){
return P_2halo_gas(ell, &z);
}else{
splin2(tab_ell-1, tab_z-1, tab_Cl_yy_2, err_Cl_yy_2, NPOINTS, NPOINTS, lnell, lnz, &pow_index);
return pow(10.0,pow_index);
}
}

void stack_yl_int(char *file){
fprintf(stderr, "reading %s...\n", file);
FILE *fp;
int i,j;

fp = fopen(file, "r");
if(fp == NULL){
fprintf(stderr, "not found %s\n", file);
}

fscanf(fp, "%d\n",&NP_YL);
tab_yl_int=dmatrix(1,NP_YL,1,NP_YL);
err_yl_int=dmatrix(1,NP_YL,1,NP_YL);

for(i=0;i<NP_YL;i++){
for(j=0;j<NP_YL;j++){
fscanf(fp, "%lf %lf %lf\n", &tab_l_ls[i], &tab_xout[j], &tab_yl_int[i+1][j+1]);
}
}

for(i=0;i<NP_YL;i++){
for(j=0;j<NP_YL;j++){
tab_yl_int[i+1][j+1] = log10(tab_yl_int[i+1][j+1]);
}
}

splie2(tab_l_ls-1, tab_xout-1,tab_yl_int,NP_YL,NP_YL,err_yl_int);
printf("set spline yl_integral...\n");
}

double y_l_integral_fast(double l_ls, double xout){
double lnell = log10(l_ls); 
double pow_index;

if(xout > tab_xout[NP_YL-1] || xout < tab_xout[0]){
printf("you can not spline y_l_integral!\n");
printf("l_ls = %e, xout= %e\n", l_ls, xout);
exit(1);
}
if(l_ls > pow(10., tab_l_ls[NP_YL-1])){
double log_lls1 = tab_l_ls[NP_YL-1];
double log_lls2 = tab_l_ls[NP_YL-2];

double res1, res2;
splin2(tab_l_ls-1, tab_xout-1, tab_yl_int, err_yl_int, NP_YL, NP_YL, log_lls1, xout, &res1);
splin2(tab_l_ls-1, tab_xout-1, tab_yl_int, err_yl_int, NP_YL, NP_YL, log_lls2, xout, &res2);

pow_index = (res1-res2)/(log_lls1-log_lls2)*(lnell-log_lls2)+res2;
return pow(10., pow_index);

}
if(l_ls < pow(10., tab_l_ls[0])){
double log_lls1 = tab_l_ls[0];
double log_lls2 = tab_l_ls[1];

double res1, res2;
splin2(tab_l_ls-1, tab_xout-1, tab_yl_int, err_yl_int, NP_YL, NP_YL, log_lls1, xout, &res1);
splin2(tab_l_ls-1, tab_xout-1, tab_yl_int, err_yl_int, NP_YL, NP_YL, log_lls2, xout, &res2);

pow_index = (res1-res2)/(log_lls1-log_lls2)*(lnell-log_lls2)+res2;
return pow(10., pow_index);
}

splin2(tab_l_ls-1, tab_xout-1, tab_yl_int, err_yl_int, NP_YL, NP_YL, lnell, xout, &pow_index);
return pow(10., pow_index);
}

double y_l(double ell, double z, double M){
double m500=M500(z,M)/1.2; //see 1509.05134v1.pdf  
double r500=R500_fast(z,M)/pow(1.20, 0.3333); //1509.05134v1.pdf

double Ez=H_z(z)/100.0/HubbleParam;

double rvir = r_vir(z, M);
double covd = chi_fast(z);
double ells = covd/(1+z)/r500;

double l_ls = ell/ells;
double xout = 6.0; //2.0*rvir/r500; see 1509.05134v1.pdf

double y_int = y_l_integral_fast(l_ls, xout);

y_int = y_int * pow(HubbleParam, -3/2);

double eps_p, A;//cf.arXiv:1106.3208 Fig 4
//TBO2
//eps_p=0.80;
//A = 1.26;

//Battaglia
//eps_p=0.08;
//A = 0.99;

//Arnaud
eps_p =0.;
A = 1.;

double res;
double index = 0.787;
res = sqrt(A)*1.65/0.7/0.7*pow(m500/3e14/0.7,index)*y_int*pow(Ez,8./3.-eps_p)*HubbleParam*HubbleParam; //eV/cc
if(isnan(res)){printf("nan Press profile at z=%e M=%e\n",z,M);}

double Mpc2cm = 3.0856*1e18*1e6;

res = res * sigma_T/m_elect * Mpc2cm; // Mpc^-1

return 4*PI*(r500/HubbleParam)/ells/ells*res; //dimensionless

}
*/

//mass-mass
double U_nfw(double M, void *params){
	double res;
	struct nfw_params *p = (struct nfw_params *)params;
	
	double k = p->wave_num;
	double cnfw = p->con;
	double rs = p->rad;
	double rhos = p->rho;
	
	double Si,Co;
	double dummy = 0,abserr;
	double x_lo = k*rs;
	double x_hi = (1.+cnfw)*k*rs;
	
	Si = Sin_Integral_Si(x_hi)-Sin_Integral_Si(x_lo);
	Co = Cos_Integral_Ci(x_hi)-Cos_Integral_Ci(x_lo);
	
	//res = integral/integral_norm;
	//k = p->wave_num;
	//if( k <= 1.e-3){res = 1.;}
	
	double integral = sin(x_lo)*Si- sin(cnfw*x_lo)/x_hi + cos(x_lo)*Co;
	res = integral/(log(1.+cnfw)-cnfw/(1.+cnfw));
	
	return res;
}

double P_nonlinear(double z, double k){
	
	double ksig ,n_eff, C_halo;
	double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
	double param[3];
	set_halofit_param(z, param);
	ksig = param[0];
	n_eff = param[1];
	C_halo = param[2];
	double a_n,b_n,c_n,alpha_n,gamma_n,beta_n,nu_n,mu_n;
	if(OPT_fit==1){
		a_n = pow(10., 1.4861 + 1.8369*n_eff + 1.6762*n_eff*n_eff + 0.7940*n_eff*n_eff*n_eff + 0.1670*n_eff*n_eff*n_eff*n_eff -0.6206*C_halo);
		b_n = pow(10.,0.9463 + 0.9466*n_eff + 0.3084*n_eff*n_eff - 0.9400*C_halo);
		c_n = pow(10.,-0.2807 + 0.6669*n_eff + 0.3214*n_eff*n_eff - 0.0793*C_halo);
		gamma_n = 0.8649 + 0.2989*n_eff + 0.1631*C_halo;
		alpha_n = 1.3884 + 0.3700*n_eff - 0.1452*n_eff*n_eff;
		beta_n = 0.8291 + 0.9854*n_eff + 0.3401*n_eff*n_eff;
		mu_n = pow(10.,-3.5442 + 0.1908*n_eff);
		nu_n = pow(10., 0.9589 + 1.2857*n_eff);
	}
	if(OPT_fit==2){
		double om_w = OmegaLambda*pow(1+z,3*w)/(Omega+OmegaLambda*pow(1+z,3*w)); 
		a_n = pow(10.,1.5222 + 2.8553*n_eff + 2.3706*n_eff*n_eff + 0.9903*n_eff*n_eff*n_eff + 0.2250*n_eff*n_eff*n_eff*n_eff -0.6038*C_halo +0.1749*om_w*(1.+w));
		b_n = pow(10.,-0.5642 + 0.5864*n_eff + 0.5716*n_eff*n_eff - 1.5474*C_halo +0.2279*om_w*(1.+w));
		c_n = pow(10.,0.3698 + 2.0404*n_eff + 0.8161*n_eff*n_eff + 0.5869*C_halo);
		gamma_n = 0.1971 - 0.0843*n_eff + 0.8460*C_halo;
		alpha_n = fabs(6.0835 + 1.3373*n_eff - 0.1959*n_eff*n_eff - 5.5274*C_halo);
		beta_n = 2.0379 - 0.7354*n_eff + 0.3157*n_eff*n_eff + 1.2490*n_eff*n_eff*n_eff + 0.3980*n_eff*n_eff*n_eff*n_eff - 0.1682*C_halo;
		mu_n = 0;
		nu_n = pow(10., 5.2105 + 3.6902*n_eff);
	}
	//printf("%e %e %e %e %e %e %e\n",a_n,b_n,c_n,alpha_n,gamma_n,beta_n,nu_n,mu_n);
	
	
	double delta_L,delta_H,delta_Q;
	//only support flat universe
	double omz;
	omz = Omega/(Omega+OmegaLambda*pow(1+z,3*w));
	
	
	double f1 = pow(omz,-0.0307);
	double f2 = pow(omz,-0.0585);
	double f3 = pow(omz,0.0743);
	
	double fac =4*PI*k*k*k/(2*PI)/(2*PI)/(2*PI);
	delta_L = fac*PowerSpec(k)/inv_D;
	
	double y = k/ksig;
	delta_Q = delta_L*(pow(1.+delta_L,beta_n)/(1.+alpha_n*delta_L))*exp(-y/4-y*y/8);
	
	delta_H = a_n*pow(y,3*f1)/(1+b_n*pow(y,f2)+pow(c_n*f3*y,3-gamma_n));
	delta_H = delta_H/(1+mu_n/y+nu_n/y/y);
	
	return (delta_Q+delta_H)/fac;
}

void set_halofit_param(double z, double *param){
	int i;
	double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
	thre = inv_D;
	//double init = pow(10.,-1.065+4.332e-1*(1.+z)-2.516e-2*pow(1.+z,2)+9.069e-4*pow(1.+z,3));
	//double R0=solver(z);
	double R0=solver(z);
	//check the solver
	//double dummy=0;
	//printf("check_solver %e\n",sigma2_gauss(log(R0),&dummy));
	
	param[0] = 1./R0; //k_sig
	param[1] = neff(R0);
	param[2] = C_halofit(R0);
	
}

double solver(double z){
	int status; 
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.;
	double x_lo = -5.0, x_hi = 5.0; 
	gsl_function F; 
	double params = z;
	F.function = &sigma2_gauss; 
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
		// printf("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r,  x_hi-x_lo);
	} while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	
	double R0 = exp(r);
	
	//printf("Rsig = %e\n",R0);
	return R0;
}

double get_delta_k(double k){
	double lk;
	double pow_index;
	double fac = k*k*k*4*PI/(2*PI)/(2*PI)/(2*PI);
	
	return fac*PowerSpec(k);
	
}

double sigma2_gauss_int(double lnk, void *params){
	double res;
	double R = *(double *)params;
	
	double k = exp(lnk);
	
	res = get_delta_k(k)*exp(-k*k*R*R);
	return res;
}

double sigma2_gauss(double lnR, void *params){
	double result,abserr;
	double params_int=exp(lnR);
	double z = *(double *)params;
	double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
	size_t neval;
	gsl_function F;
	
	//gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);
	
	//F.function=&sigma2_gauss_int;
	//F.params=&params_int;
	//gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
	//gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
	//gsl_integration_workspace_free(wgsl);
	
	/*int Nint = 5000;
	
	double dlnk = 10.0/Nint;
	result = 0;
	for(int i=0;i<Nint;i++){
	result += dlnk*sigma2_gauss_int(i*dlnk-5.0,&params_int);
	}*/
	//res = qromb(integral_Pkappa,log(calc_z),log(zLSS));
	
	return sigma2_gauss_fast(params_int)-inv_D;
}

double dsig_dR_int(double lnk, void *params){
	double res;
	double R = *(double *)params;
	
	double k = exp(lnk);
	
	res = get_delta_k(k)*exp(-k*k*R*R)*(-2*k*k*R);
	return res;
}
double dsig_dR(double R){
	double result,abserr;
	double params=R;
	size_t neval;
	gsl_function F;
	
	//gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);
	
	//F.function=&dsig_dR_int;
	//F.params=&params;
	//gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
	//gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
	//gsl_integration_workspace_free(wgsl);
	
	int Nint = 5000;
	
	double dlnk = 20.0/Nint;
	result = 0;
	for(int i=0;i<Nint;i++){
		result += dlnk*dsig_dR_int(i*dlnk-10.0,&params);
	}
	//res = qromb(integral_Pkappa,log(calc_z),log(zLSS));
	
	return result;
}
double d2sig_dR2_int(double lnk, void *params){
	double res;
	double R = *(double *)params;
	
	double k = exp(lnk);
	double y = k*R;
	
	res = get_delta_k(k)*exp(-y*y)*(y*y-y*y*y*y);
	
	return res;
}
double d2sig_dR2(double R){
	double result,abserr;
	double params=R;
	size_t neval;
	gsl_function F;
	
	//gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);
	
	//F.function=&d2sig_dR2_int;
	//F.params=&params;
	//gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
	//gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
	//gsl_integration_workspace_free(wgsl);
	
	int Nint = 5000;
	
	double dlnk = 20.0/Nint;
	result = 0;
	for(int i=0;i<Nint;i++){
		result += dlnk*d2sig_dR2_int(i*dlnk-10.0,&params);
	}
	//res = qromb(integral_Pkappa,log(calc_z),log(zLSS));
	//if(result < 0){printf("R=%e d^2sig/dR^2 = %e\n",R,result);}
	
	return result;
}
double neff(double R){
	double sig2 = thre;
	double res = -R*dsig_dR_fast(R)/sig2;
	
	return res -3.0;
}

double C_halofit(double R){
	double sig2 = thre;
	double n_eff = neff(R);
	
	double res = (3.+n_eff)*(3+n_eff) +4./thre*d2sig_dR2_fast(R); //cf .smith et al
	
	return res;
	
}

void stack_Pk_nonl_data_and_spline(){
	double logR = (5.-(-5.))/NPOINTS;
	int i;
	
	for(i=0;i<NPOINTS;i++){
		tab_R[i] = (double)(i)*logR -5.;
		tab_dsdr[i] = dsig_dR(pow(10.,tab_R[i]));
		tab_ds2dr2[i] = d2sig_dR2(pow(10.,tab_R[i]));
		//printf("%e %e %e\n",tab_R[i],tab_dsdr[i],tab_ds2dr2[i]);
		tab_dsdr[i] = log10(-tab_dsdr[i]);
		tab_ds2dr2[i] = log10(10+tab_ds2dr2[i]);
		//printf("%e %e %e\n",tab_R[i],tab_dsdr[i],tab_ds2dr2[i]);
	}
	
	int Nint = 5000;
	
	double dlnk = 20.0/Nint;
	double params_int;
	for(i=0;i<NPOINTS;i++){
		tab_sig2[i]=0.;
		params_int = pow(10.,tab_R[i]);
		for(int j=0;j<Nint;j++){
			tab_sig2[i] += dlnk*sigma2_gauss_int(j*dlnk-10.0,&params_int);
		}
		tab_sig2[i] = log10(tab_sig2[i]);
	}
	
	double yp1 = 1.e31;
	double ypn = 1.e31;
	
	spline(tab_R-1, tab_sig2-1, NPOINTS, yp1, ypn, err_sig2-1);
	printf("spline sigma2 ... \n");
	spline(tab_R-1, tab_dsdr-1, NPOINTS, yp1, ypn, err_dsdr-1);
	printf("spline dsig/dr ... \n");
	spline(tab_R-1, tab_ds2dr2-1, NPOINTS, yp1, ypn, err_ds2dr2-1);
	printf("spline d2sig/dr2 ... \n");
}

double dsig_dR_fast(double R){
	double lR;
	double pow_index;
	
	lR=log10(R);
	
	if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
		return dsig_dR(R);
	}else{
		splint(tab_R-1, tab_dsdr-1, err_dsdr-1, NPOINTS, lR, &pow_index);
		return -pow(10.,pow_index);
	}
}
double d2sig_dR2_fast(double R){
	double lR;
	double pow_index;
	
	lR=log10(R);
	
	if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
		return d2sig_dR2(R);
	}else{
		splint(tab_R-1, tab_ds2dr2-1, err_ds2dr2-1, NPOINTS, lR, &pow_index);
		return pow(10.,pow_index)-10.;
	}
}

double sigma2_gauss_fast(double R){
	double lR;
	double pow_index;
	
	lR=log10(R);
	
	if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
		int Nint = 5000;
		double dlnk = 20.0/Nint;
		double result=0.;
		
		for(int i=0;i<Nint;i++){
			result += dlnk*sigma2_gauss_int(i*dlnk-10.0,&R);
		}
		return result;
	}else{
		splint(tab_R-1, tab_sig2-1, err_sig2-1, NPOINTS, lR, &pow_index);
		return pow(10.,pow_index);
	}
}

/*
double mean_density_int(double logm, void *params){
double res;
double rm = pow(10,logm);

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(calc_z, rm);
logm = log10(m180m);
}

double mf = dndlogm_fast(logm, calc_z);

res = pow(10., logm)*mf;

return res;
}

double mean_density(){
double res,abserr,dummy=0;
double x_lo = 9.0;
double x_hi = log10(Mhalo_max);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(mean_density_int(x_hi, &dummy)+mean_density_int(x_lo, &dummy));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*mean_density_int(x_lo+(2*j-1)*h[i-1], &dummy);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

//printf("z=%3.2f\n", calc_z);
//printf("mean density of DM = %e [(M_sun/h)/(Mpc/h)^3]\n",res);
//printf("true density of DM = %e [(M_sun/h)/(Mpc/h)^3]\n",2.7754e11*Omega*pow(1.+calc_z,3));
return res;
}

double A00(double logm,double k){//cf.astro-ph/1009.0597v2
double res;
double rm = pow(10,logm);

double q = pow(3*rm/4/PI/rho_mean,0.3333333333);
double kq = k*q;

res = 3*(sin(kq)/kq/kq/kq - cos(kq)/kq/kq);

return res;
}

double shot_noise_int(double logm, void *params){
double res;
double calc_k = *(double *)params;
double rm = pow(10,logm);

double mf = dndlogm_fast(logm, calc_z);
double ceff = A00(logm, calc_k);

res = pow(10., logm)*pow(10., logm)*mf*ceff*ceff;
return res;
}

double shot_noise(double k){
double res,abserr,dummy=0;
double x_lo = 9.0;
double x_hi = log10(Mhalo_max);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(shot_noise_int(x_hi, &k)+shot_noise_int(x_lo, &k));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*shot_noise_int(x_lo+(2*j-1)*h[i-1], &k);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];
return res;
}

double onehalo_int_mm(double logm, void *params){
double rm = pow(10.0,logm);

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(calc_z, rm);
logm = log10(m180m);
}

double mf = dndlogm_fast(logm, calc_z);

double calc_k = *(double *)params;
double cnfw, rs, rhos;

cnfw = c_nfw(calc_z, rm);
rs =   r_s(calc_z, rm);
rhos = rho_s(calc_z, rm);

struct nfw_params p = {calc_k,cnfw,rs,rhos};

double u = U_nfw(rm, &p);
double res = mf*pow(10., logm)*pow(10., logm)*u*u;

return res;

}
double twohalo_int_mm(double logm, void *params){
double rm = pow(10,logm);

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(calc_z, rm);
logm = log10(m180m);
}

double mf = dndlogm_fast(logm, calc_z);

double calc_k = *(double *)params;
double cnfw, rs, rhos;
cnfw = c_nfw(calc_z,rm);
rs   = r_s(calc_z,  rm);
rhos = rho_s(calc_z,rm);

struct nfw_params p = {calc_k,cnfw,rs,rhos};

double u = U_nfw(rm, &p);
double b = halo_bias_fast(logm, calc_z);

double res = mf*pow(10., logm)*b*u;

return res;
}

double P_1halo_mm(double k){
double res,abserr;
double x_lo = 9.0;
double x_hi = log10(Mhalo_max);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(onehalo_int_mm(x_hi, &k)+onehalo_int_mm(x_lo, &k));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*onehalo_int_mm(x_lo+(2*j-1)*h[i-1], &k);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}

double P_2halo_mm(double k){
double res,abserr;
double x_lo = 9.0;
double x_hi = log10(Mhalo_max);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(twohalo_int_mm(x_hi, &k)+twohalo_int_mm(x_lo, &k));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*twohalo_int_mm(x_lo+(2*j-1)*h[i-1], &k);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

double gfac=(growth(1./(1+calc_z))/growth(1.0));
return res*res*PowerSpec(k)*gfac*gfac;
}

void stack_Pks(int opt){
int i,j;
double dlogk = (5.0-(-4.0))/NPOINTS;

tab_pdd1=dmatrix(1,NPOINTS,1,NPOINTS);
err_pdd1=dmatrix(1,NPOINTS,1,NPOINTS);
tab_pdd2=dmatrix(1,NPOINTS,1,NPOINTS);
err_pdd2=dmatrix(1,NPOINTS,1,NPOINTS);

for(i=0;i<NPOINTS;i++){tab_k[i] = -4.0+(double)(i)*dlogk;}

for(i=0;i<NPOINTS;i++){
fprintf(stderr, ".");
for(j=0;j<NPOINTS;j++){
double calc_k = pow(10, tab_k[i]);
calc_z = pow(10, tab_z[j]);
//double mean1 = 2.7754e11*Omega*(1.+calc_z)*(1.+calc_z)*(1.+calc_z);
double mean2 = mean_density();
//rho_mean = mean;
tab_pdd1[i+1][j+1] = P_1halo_mm(calc_k)/pow(mean2,2);
tab_pdd2[i+1][j+1] = P_2halo_mm(calc_k)/pow(mean2,2);

if(tab_pdd1[i+1][j+1] < 1e-30 || tab_pdd2[i+1][j+1] < 1e-30){
fprintf(stdout, "%e %e %e %e\n", 
calc_k, calc_z, tab_pdd1[i+1][j+1], tab_pdd2[i+1][j+1]);
exit(1);
}

tab_pdd1[i+1][j+1] = log10(tab_pdd1[i+1][j+1]);
tab_pdd2[i+1][j+1] = log10(tab_pdd2[i+1][j+1]);
}
}

splie2(tab_k-1, tab_z-1, tab_pdd1, NPOINTS, NPOINTS, err_pdd1);
splie2(tab_k-1, tab_z-1, tab_pdd2, NPOINTS, NPOINTS, err_pdd2);

cout << "calc Pk's..." << endl;
}

double P_1halo_mm_fast(double k, double z){
double lk = log10(k);
double lz = log10(z);
double pow_index;

if(k < pow(10., tab_k[0]) || k > pow(10., tab_k[NPOINTS-1]))return 0;
if(z < pow(10., tab_z[0]) || z > pow(10., tab_z[NPOINTS-1]))return 0;

splin2(tab_k-1, tab_z-1, tab_pdd1, err_pdd1, NPOINTS, NPOINTS, lk, lz, &pow_index);
//printf("OK.spline dn/dlogm\n");
if(isnan(pow_index)){
printf("fail to spline P_dd1 at k=%e z=%e\n",k,z);
exit(1);
}
return pow(10.0,pow_index);
}
double P_2halo_mm_fast(double k, double z){
double lk = log10(k);
double lz = log10(z);
double pow_index;

if(k < pow(10., tab_k[0]) || k > pow(10., tab_k[NPOINTS-1]))return 0;
if(z < pow(10., tab_z[0]) || z > pow(10., tab_z[NPOINTS-1]))return 0;

splin2(tab_k-1, tab_z-1, tab_pdd2, err_pdd2, NPOINTS, NPOINTS, lk, lz, &pow_index);
//printf("OK.spline dn/dlogm\n");
if(isnan(pow_index)){
printf("fail to spline P_dd2 at k=%e z=%e\n",k,z);
exit(1);
}
return pow(10.0,pow_index);
}
*/

//lensing staff
double da_comb(double zl){
	double result,abserr;
	size_t neval;
	gsl_function F;
	
	if(OPT_survey == 0){
		double dl,ds,dls;
		dl = chi(zl);
		ds = chi(z_source);
		if(z_source >= zl)dls = ds-dl;
		if(z_source < zl)dls = 0.0;
		result = dls/ds;
		goto SKIP;
	}
	
	F.function=&da_comb_int;
	F.params=&zl;
	if(zl > 1e-5){
		//gsl_integration_qng(&F,log(zl),+5,0,1e-6, &result, &abserr,&neval);
		
		gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);
		gsl_integration_qag(&F,log(zl),log(10.), 0, 1e-4, 1000, 6,wgsl, &result, &abserr);
		gsl_integration_workspace_free(wgsl);
		
	}else{
		result=0;
	}
	
	result = result/gal_norm;
	
  SKIP:
	//printf("%e %e\n",zl,result);
	return result;
}

double da_comb_int(double lnz, void *param){
	double res;
	double zl = *(double *)param;
	double z = exp(lnz);
	double dummy=0;
	double dl,ds,dls;
	
	dl  = chi(zl);
	ds  = chi(z);
	dls = ds-dl;
	
	res = z*gal_dist(z)*dls/ds;
	
	return res;
}

double gal_dist(double z){
	double res;
	double z_median,z0;
	double a_z,b_z,c_z,d_z;
	z_median = z_source;
	int iz;
	
	switch(OPT_survey)
	{
	  case 1://cf. 1010.0744
		z0 = z_median/3.;
		res = z*z*exp(-(z/z0))/(2*z0*z0*z0);
		return res;
		break;
		
	  case 2://cf. 1204.4530
		z0=0.71*z_median;
		res=3*z*z*exp(-pow(z/z0,1.5))/(2*z0*z0*z0);
		return res;
		break;
		
	  case 3://cf. 0712.0884
		a_z=0.6;
		b_z=7.0;
		c_z=0.6;
		if(z<1e-5){
			res=0.;
		}else{
			res=(pow(z,a_z)+pow(z,a_z*b_z))/(pow(z,b_z)+c_z);
		}
		return res;
		break;
	  case 4://cf.http://xxx.lanl.gov/pdf/1303.1806v2.pdf
		a_z = 1.50;
		b_z = 0.32;
		c_z = 0.20;
		d_z = 0.46;
		res = a_z * exp(-(z-0.7)*(z-0.7)/b_z/b_z)+c_z*exp(-(z-1.2)*(z-1.2)/d_z/d_z);
		return res;
		break;
	  default:
		fprintf(stdout,"gal_dist Not supported\n");  
	}
}

double gal_dist_int(double lnz, void *param){
	double res;
	double z=exp(lnz);
	
	res = gal_dist(z)*z;
	return res;
}

double gal_dist_norm(double zmax){
	double result,abserr;
	size_t neval;
	
	double param=0;
	
	gsl_function F;
	
	F.function=&gal_dist_int;
	F.params=&param;
	//gsl_integration_qng(&F, -5., log(zmax), 0, 1e-6, &result, &abserr, &neval);
	gsl_integration_workspace *w_gsl= gsl_integration_workspace_alloc(1000);
	gsl_integration_qag(&F, -5., log(zmax), 0, 1e-6, 1000, 6, w_gsl, &result, &abserr);
	gsl_integration_workspace_free(w_gsl);
	
	return result;
}

double Window_kappa(double z){
	double res;
	res = 1.5*(100/C)*(100/C)*Omega*(1.+z)*chi_fast(z)*da_comb(z);
	return res;
}

double lens_kernel(double z){//only support flat universe
	double res;
	double fac = 1.5*(100/C)*(100/C)*Omega;
	//res = fac*(chi(calc_z)-chi(z))/chi(calc_z)*(growth(1./(1+z))/growth(1.0))*(1+z);
	
	res = fac*da_comb(z)*(1.+z);
	res = res*res;
	
	return C*res*HubbleParam/H_z(z);
}

double integral_Pkappa(double lnz, int opt){
	double res;
	double z = exp(lnz);
	double k = multipole/chi_fast(z);
	double gfac = growth(1./(1+z))/growth(1.0);
	switch(opt)
	{
	  case 1:
		res=z*lens_kernel(z)*PowerSpec(k)*gfac*gfac;
		break;
	  case 2:
		if(z <= 3){res = z*lens_kernel(z)*P_nonlinear(z,k);}
		if(z > 3){res =  z*lens_kernel(z)*PowerSpec(k)*gfac*gfac;}
		break;
	  default:
		fprintf(stdout,"llcl Not supported\n");  
	}
	
	//if(res > 1e-30)printf("z=%e %e\n",z,res);
	return res;
	
	/*
	double res;
	double z = exp(lnz);
	double k = multipole/chi_fast(z);
	
	switch(opt)
	{
	case 1:
	res=z*lens_kernel(z)*P_1halo_mm_fast(k, z);
	break;
	case 2:
	res=z*lens_kernel(z)*P_2halo_mm_fast(k, z);
	break;
	default:
	fprintf(stdout,"llcl Not supported\n");  
	}
	*/
	
	//if(res > 1e-30)printf("z=%e %e\n",z,res);
	return res;
}

double Pkappa(double ell, int opt){
	double res;	
	multipole = ell;
	
	double x_lo = log(1e-3);
	double x_hi = log(6.);
	
	if(OPT_survey==0)x_hi = log(z_source);
	
	int Nint = 10;
	//Romberg
	int i,j;
	double h[Nint];
	double s[Nint][Nint];
	for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}
	
	s[0][0] = 0.5*h[0]*(integral_Pkappa(x_hi, opt)+integral_Pkappa(x_lo, opt));
	for(i=2;i<=Nint;i++){
		s[i-1][0] = s[i-2][0];
		for(j=1;j<=pow(2.,i-2);j++){
			s[i-1][0] += h[i-2]*integral_Pkappa(x_lo+(2*j-1)*h[i-1], opt);
		}
		s[i-1][0] = 0.5*s[i-1][0];
	}
	
	for(i=2;i<=Nint;i++){
		for(j=2;j<=i;j++){
			s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
		}
	}
	
	res = s[Nint-1][Nint-1];
	
	return res;
}

//cross correlation
/*
double onehalo_int_mgas(double logm, void *params){
double rm = pow(10,logm);

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(calc_z, rm);
logm = log10(m180m);
}

double mf = dndlogm_fast(logm, calc_z);

double ell = *(double *)params;
double cnfw, rs, rhos;

cnfw = c_nfw(calc_z,rm);
rs   = r_s(calc_z,rm);
rhos = rho_s(calc_z,rm);

double dl = chi_fast(calc_z);
double calc_k = ell/dl;

struct nfw_params p = {calc_k,cnfw,(1.+calc_z)*rs,rhos};

double norm_m = 4*PI*rhos*rs*rs*rs*(log(1.+cnfw)-cnfw/(1.+cnfw));
double Kappa = norm_m*U_nfw(rm, &p)*(1.+calc_z)*(1.+calc_z)/dl/dl;
double y  = y_l(ell, calc_z, rm);

double res = mf*Kappa*y;

return res;
}

double twohalo_int_mgas(double logm, void *params){
double rm = pow(10,logm);

if(OPT == 0){ // Mvir to M180m
double m180m = M180m(calc_z, rm);
logm = log10(m180m);
}

double mf = dndlogm_fast(logm, calc_z);

double ell = *(double *)params;
double cnfw, rs, rhos;

cnfw = c_nfw(calc_z,rm);
rs   = r_s(calc_z,rm);
rhos = rho_s(calc_z,rm);

double dl = chi_fast(calc_z);
double calc_k = ell/dl;

struct nfw_params p = {calc_k,cnfw,(1.+calc_z)*rs,rhos};

double norm_m = 4*PI*rhos*rs*rs*rs*(log(1.+cnfw)-cnfw/(1.+cnfw));
double Kappa = norm_m*U_nfw(rm, &p)*(1.+calc_z)*(1.+calc_z)/dl/dl;
double b = halo_bias_fast(logm, calc_z);

double res = mf*b*Kappa;

return res;
}

double P_1halo_mgas(double k){
double res,abserr;
size_t neval;
double x_lo = log10(5.0e11);
double x_hi = log10(5.0e15);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(onehalo_int_mgas(x_hi,&k)+onehalo_int_mgas(x_lo,&k));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*onehalo_int_mgas(x_lo+(2*j-1)*h[i-1],&k);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}

double P_2halo_mgas(double k){
double res1,res2, abserr;
size_t neval;
double x_lo = log10(5.0e11);
double x_hi = log10(5.0e15);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(twohalo_int_mgas(x_hi,&k)+twohalo_int_mgas(x_lo,&k));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*twohalo_int_mgas(x_lo+(2*j-1)*h[i-1],&k);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res1 = s[Nint-1][Nint-1];

double dA = chi_fast(calc_z);
double gfac=(growth(1./(1+calc_z))/growth(1.0));
struct gas_params p_int = {k,calc_z};

for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(twohalo_int_gas(x_hi,&p_int)+twohalo_int_gas(x_lo,&p_int));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*twohalo_int_gas(x_lo+(2*j-1)*h[i-1],&p_int);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res2 = s[Nint-1][Nint-1];

return res1*res2*gfac*gfac*PowerSpec(k/dA);
}

double P_1halo_mgas_fast(double ell, double z){
double lnz = log10(z);
double lnell = log10(ell); 

double pow_index;

if(lnz > tab_z[NPOINTS-1] || lnz < tab_z[0] || ell > pow(10., tab_ell[NPOINTS-1]) || ell < pow(10., tab_ell[0])){
calc_z = z;
return P_1halo_mgas(ell);
}else{
splin2(tab_ell-1, tab_z-1, tab_Cl_yk_1, err_Cl_yk_1, NPOINTS, NPOINTS, lnell, lnz, &pow_index);
return pow(10.0,pow_index);
}
}

double P_2halo_mgas_fast(double ell, double z){
double lnz = log10(z);
double lnell = log10(ell); 

double pow_index;

if(lnz > tab_z[NPOINTS-1] || lnz < tab_z[0] || ell > pow(10., tab_ell[NPOINTS-1]) || ell < pow(10., tab_ell[0])){
calc_z = z;
return P_2halo_mgas(ell);
}else{
splin2(tab_ell-1, tab_z-1, tab_Cl_yk_2, err_Cl_yk_2, NPOINTS, NPOINTS, lnell, lnz, &pow_index);
return pow(10.0,pow_index);
}
}

double Cl_1halo_yk_int(double lnx, void *params){
double res;
double ell = *(double *)params;
double z = exp(lnx); //lnx = log(z)
double covd = chi_fast(z);

double sigma_crit_inv=4*PI*G0/C/C*covd*da_comb(z)/(1.+z); //(Msun/h)^-1

if(z > 1e-4){
res = z*covd*covd*C*HubbleParam/H_z(z)*P_1halo_mgas_fast(ell, z)*sigma_crit_inv;
}else{
res=1e-30;
}
//printf("%g %g %g\n",ell,z,res);

return res;
}
double Cl_1halo_yk(double ell, double zmax){
double res,abserr;
size_t neval;
double x_lo = log(1e-3);
double x_hi = log(zmax);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(Cl_1halo_yk_int(x_hi, &ell)+Cl_1halo_yk_int(x_lo, &ell));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*Cl_1halo_yk_int(x_lo+(2*j-1)*h[i-1], &ell);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}

double Cl_2halo_yk_int(double lnx, void *params){
double res;
double ell = *(double *)params;
double z = exp(lnx); //lnx = log(z)
double covd = chi_fast(z);

double sigma_crit_inv=4*PI*G0/C/C*covd*da_comb(z)/(1.+z); //(Msun/h)^-1

if(z > 1e-4){
res = z*covd*covd*C*HubbleParam/H_z(z)*P_2halo_mgas_fast(ell, z)*sigma_crit_inv;
}else{
res=1e-30;
}
//printf("%g %g %g\n",ell,z,res);

return res;
}

double Cl_2halo_yk(double ell, double zmax){
double res,abserr;
size_t neval;
double x_lo = log(1e-3);
double x_hi = log(zmax);

int Nint = 10;
//Romberg
int i,j;
double h[Nint];
double s[Nint][Nint];
for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

s[0][0] = 0.5*h[0]*(Cl_2halo_yk_int(x_hi,&ell)+Cl_2halo_yk_int(x_lo,&ell));
for(i=2;i<=Nint;i++){
s[i-1][0] = s[i-2][0];
for(j=1;j<=pow(2.,i-2);j++){
s[i-1][0] += h[i-2]*Cl_2halo_yk_int(x_lo+(2*j-1)*h[i-1],&ell);
}
s[i-1][0] = 0.5*s[i-1][0];
}

for(i=2;i<=Nint;i++){
for(j=2;j<=i;j++){
s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
}
}

res = s[Nint-1][Nint-1];

return res;
}
*/

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

double get_eq_for_M_vir(double x, void *params){
	struct get_mvir *p = (struct get_mvir *)params;
	
	double delta = p->delta;
	double rd = p->rd; 
	double z = p->redshift;
	
	//solve 3 rho_s [ln(1+c_delta)-c_delta/(1+c_delta)] = delta * rho_m(z) * c_delta^3
	double rhos = rho_s(z, pow(10.,x));
	double rs = r_s(z, pow(10.,x));
	
	double cd = rd/rs;
	double rho_crit = 2.7754e11*(Omega*pow(1+z,3));
	
	//printf("%e %e %e %e\n",delta, delta_v(z), cd, cv);
	
	double rhs = delta*rho_crit*cd*cd*cd;
	double lfs = 3*rhos*(log(1+cd)-cd/(1+cd));
	return rhs-lfs;
}

double M_delta_to_M_vir(double z, double mdelta, double delta){
	double rd = r_delta(z, mdelta, delta);
	
	//get mvir
	int status; 
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 1.0;
	double x_lo = 5.0, x_hi = 19.0;
	gsl_function F;
	struct get_mvir params={delta, rd, z};
	
	//printf("%e %e\n",get_eq_for_M_vir(x_lo, &params), get_eq_for_M_vir(x_hi, &params));
	
	F.function = &get_eq_for_M_vir;
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
		
	} while (status == GSL_CONTINUE /*&& iter < max_iter*/);
	
	gsl_root_fsolver_free(s);
	
	return pow(10.,r);
}

double M_delta_to_M_vir_fast(double z, double mdelta, double delta){
	double lz = log10(z);
	double logm = log10(mdelta);
	
	double pow_index;
	
	if(lz > tab_z[NPOINTS-1] || lz < tab_z[0] || logm > tab_m[NPOINTS-1] || logm < tab_m[0]){
		return M_delta_to_M_vir(z, mdelta, delta);
	}else{
		splin2(tab_z-1, tab_m-1, tab_mdeltatomvir, err_mdeltatomvir, NPOINTS, NPOINTS, lz, logm, &pow_index);
		//printf("OK.spline dn/dlogm\n");
		if(isnan(pow_index))printf("fail to spline at m=%e z=%e\n",pow(10.,logm),z);
		return pow(10.0,pow_index);
	}
	
}

double get_eq_for_M_delta(double x, void *params){
	struct get_mvir *p = (struct get_mvir *)params;
	
	double delta = p->delta;
	double mvir = p->rd; //this is mvir, here
	double z = p->redshift;
	
	//solve delta/delta_vir f(1/c) = f(1/c rvir/rdelta) where f(x)=x^3 [ln(1+1/x)-1/(1+x)]
	double cv = c_nfw(z, mvir);
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
	double rs = r_s(z, mvir);
	double rd = cd*rs;
	double rho_crit=2.7754e11*(Omega*pow(1+z,3)/*+OmegaLambda*/);
	double mdelta = 4*PI/3*delta*rho_crit*rd*rd*rd;
	
	return mdelta;
}

//integral with input tables
double tab_spline_and_integral(int Nbin, double *xlist, double *ylist, double *zlist){
	int i,j;
	double x_hi = xlist[Nbin-1];
	double x_lo = xlist[0];
	int Nint = 10;
	//Romberg
	double h[Nint];
	double s[Nint][Nint];
	for(int ii=1;ii<=Nint;ii++){h[ii-1] = (x_hi-x_lo)/pow(2.,ii-1);}
	
	s[0][0] = 0.5*h[0]*(ylist[Nbin-1]+ylist[0]);
	for(int ii=2;ii<=Nint;ii++){
		s[ii-1][0] = s[ii-2][0];
		for(j=1;j<=pow(2.,ii-2);j++){
			double res, logx = x_lo+(2*j-1)*h[ii-1];
			splint(xlist-1, ylist-1, zlist-1, Nbin, logx, &res);
			s[ii-1][0] += h[ii-2]*res;
		}
		s[ii-1][0] = 0.5*s[ii-1][0];
	}
	
	for(int ii=2;ii<=Nint;ii++){
		for(j=2;j<=ii;j++){
			s[ii-1][j-1] = s[ii-1][j-2]+(s[ii-1][j-2]-s[ii-2][j-2])/(pow(4.,j-1)-1);
		}
	}
	
	return s[Nint-1][Nint-1];  
}

double chi2Redshift(double x,double z_in){
	//x,x(z_in) [Mpc/h]  dx [Mpc]  
	int i;
	int nint = 5000;
	double dx=(x-chi_fast(z_in))/HubbleParam/nint;
	double z=z_in,dz;
	
	for(i=0;i<nint;i++){
		dz = H_z(z)/C*dx;
		z += dz;
	}
	return z;
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
	double b_ling=GrowthFactor(1.,1./(1.+z));
	Const_MF = b_ling*b_ling;
	double sig = sigma_m(rmass, &r);
	
	return 1.686/sig;
}

double nu_M_fast(double z, double M){
	double logm = log10(M);
	double lz = log10(z);
	double pow_index;
	
	splin2(tab_z-1, tab_m-1, tab_nuM, err_nuM, NPOINTS, NPOINTS, lz, logm, &pow_index);
	//printf("OK.spline dn/dlogm\n");
	if(isnan(pow_index)){
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
	if(isnan(pow_index)){
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
