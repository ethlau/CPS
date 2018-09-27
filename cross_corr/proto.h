#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

double chi2Redshift(double x,double z_in);

double PowerSpec(double kmag);
double GrowthFactor(double astart, double aend);
double growth(double a);
inline double growth_int(double);
double qromb(double (*func)(double), double a, double b);
double sigma2_int(double k);
double xi_int(double k);
inline double TopHatSigma2(double R);
inline double xi(double R);
inline double PowerSpec_Efstathiou(double k);
inline double PowerSpec_EH(double k);
inline double PowerSpec_EH_neutrino(double k);
inline double PowerSpec_BBKS(double k);
inline double PowerSpec_CMBFAST(double k);
inline double PowerSpec_2D(double k);


int initialize_powerspectrum(int Spectrum);
int set_units();

double   tk_neutrino(double k);
double   growthD1(double a);
double   tk_eh(double k);
double transfunc_cmbfast(double k);
double transfunc_WDM(double k);
inline double integrand_minus1(double k);
inline double integrand_P(double k);
inline double integrand_0(double k);
inline double integrand_1(double k);
inline double integrand_2(double k);
inline double integrand_3(double k);
inline double WTopHat(double x);
inline double WGaussian(double x);
inline double F_Omega(double a);
inline double Hubble_a(double a);
double Window(double x);

inline double delta_c_func();
inline double Press_Schechter(double sigma);
inline double Sheth_Tormen(double sigma);
inline double efn(double x, double);
inline double var3(double x, double);
inline double dweight(double x, double);
inline double weight(double x, double);
double dlogdsigma(double mass, double, double);
double sigma_m(double m, double *rsphere_return);
double fract_mass(double sig);
double sigdsigdr(double);
double unnsigma(double);
inline double evar2(double x, double);
inline double var2(double x, double);
int output_scatter(char *name, double *target1, double *target2, int &nsize);
double sigmam_ks(double mass, double omega0, double hubble,double gamma, double sigma8);
double dndlogm(double logm);
double g_multi(double logm);
double sigmam_linear_eh(double mass);
double dsdm_eh(double mass);
double n_ps(double mass, double zred);
double n_st(double mass, double zred);
double n_jk(double mass, double zred);
inline double n_RSI(double k);


void readCMB_and_do_spline();
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);


double Omega_de(double a);
double coeff1(double a);
double coeff2(double a);//u"+coeff1(a)u'+coeff2(a)u=0, u=D/a, '=d/dlna
double RungeKutta(double a_in,double a); //calculate linear density growth eq.
void growth_spline();
double growth_for_any_w(double a);

void stack_table_and_spline();
void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y);

double chi_fast(double z);
double dndlogm_fast(double logm, double z);
double halo_bias_fast(double logm, double z);
double R500_fast(double z, double M);
double rvir_fast(double z, double M);

double conv_Press_fast(double M, void *params);

double halo_bias(double logm);
double H_z(double z);
double H_z_1(double z,void *params);
double integral_chi(double lnz,void *params);
double chi(double z);
double rho_nfw(double r, double z, double M);
double rho_s(double z,double M);
double r_s(double z,double M);
double c_nfw(double z,double M);
double delta_c(double z,double M);
double delta_v(double z);
double r_vir(double z,double M);

double M500(double z, double M);
double R500(double z, double M);
double R500_v2(double z,double M);

double calc_R500_func(double x, void *params);
double Rstar(double z, double M);
double calc_Rstar_func(double x, void *params);
double Psi_nfw(double r,double z,double M);
double sigma_v(double r,double z,double M);
double Li2(double x);
double totE_nfw(double r,double z, double M);

double R180m(double z, double M);
double R180m_fast(double z, double M);
double M180m(double z, double M);

double R200(double z, double M);
double R200_fast(double z, double M);
double M200(double z, double M);

//compton-y statistics
double onehalo_int_gas(double logm, void *params);
double twohalo_int_gas(double logm, void *params);

double P_1halo_gas(double k, void *params);
double P_2halo_gas(double k, void *params);

double Cl_1halo_yy_int(double lnx, void *params);
double Cl_1halo_yy(double ell, double zmax);

double Cl_2halo_yy_int(double lnx, void *params);
double Cl_2halo_yy(double ell, double zmax);

void stack_conv_Press();
double P_1halo_gas_fast(double ell, double z);
double P_2halo_gas_fast(double ell, double z);

void stack_yl_int(char *file);
double y_l_integral_int_x(double x, double l_ls);
double y_l_integral(double l_ls, double xout);
double y_l_integral_fast(double l_ls, double xout);
double y_l(double ell, double z, double M);

double P200_REXCESS_int(double x, double z, double M);
double P200_REXCESS(double z, double M);

double P200_B12_int(double x, double z, double M);
double P200_B12(double z, double M); //http://arxiv.org/pdf/1412.5593.pdf

double U_nfw(double M, void *params);

//1halo term
//                         Internally Defined Routines                        //
double      Sin_Integral_Si( double x );
double      Entire_Cos_Integral_Cin( double x );
double      Cos_Integral_Ci( double x );
void        Sin_Cos_Integrals_Si_Ci( double x, double *Si, double *Ci );
long double xSin_Integral_Si( long double x );
long double xEntire_Cos_Integral_Cin( long double x );
long double xCos_Integral_Ci( long double x );
void        xSin_Cos_Integrals_Si_Ci( long double x, long double *Si,
			long double *Ci );

static long double Asymptotic_Series_Ci( long double x );

//                         Externally Defined Routines                        //
extern void xAuxiliary_Sin_Cos_Integrals_fi_gi(long double x, long double *fi, 
			long double *gi);
extern long double xPower_Series_Si( long double x );
extern long double xPower_Series_Cin( long double x );

double mean_density_int(double logm, void *params);
double mean_density();
double A00(double logm,double k);//cf.astro-ph/1009.0597v2
double shot_noise_int(double logm, void *params);
double shot_noise(double k);
double onehalo_int_mm(double logm, void *params);
double twohalo_int_mm(double logm, void *params);
double P_1halo_mm(double k);
double P_2halo_mm(double k);
void stack_Pks(int opt);
double P_1halo_mm_fast(double k, double z);
double P_2halo_mm_fast(double k, double z);

double P_nonlinear(double z, double k);
void set_halofit_param(double z, double *param);
double solver(double z);
double get_delta_k(double k);
double sigma2_gauss_int(double lnk, void *params);
double sigma2_gauss(double lnR, void *params);
double dsig_dR_int(double lnk, void *params);
double dsig_dR(double R);
double d2sig_dR2_int(double lnk, void *params);
double d2sig_dR2(double R);
double neff(double R);
double C_halofit(double R);
void stack_Pk_nonl_data_and_spline();
double dsig_dR_fast(double R);
double d2sig_dR2_fast(double R);
double sigma2_gauss_fast(double lnR);


//lensing staff
double da_comb(double z);
double da_comb_int(double lnz, void *param);
double gal_dist(double z);
double gal_dist_int(double lnz,void *param);
double gal_dist_norm(double zmax);
double Window_kappa(double z);
double lens_kernel(double z);
double integral_Pkappa(double lnz, int opt);
double Pkappa(double ell, int opt);

//cross correlation
double onehalo_int_mgas(double logm, void *params);
double twohalo_int_mgas(double logm, void *params);

double P_1halo_mgas(double k);
double P_2halo_mgas(double k);

double P_1halo_mgas_fast(double k, double z);
double P_2halo_mgas_fast(double k, double z);

double Cl_1halo_yk_int(double lnx, void *params);
double Cl_1halo_yk(double ell, double zmax);

double Cl_2halo_yk_int(double lnx, void *params);
double Cl_2halo_yk(double ell, double zmax);

// Mvir <-> Mdelta
double r_delta(double z, double Mass, double delta);
double get_eq_for_M_vir(double x, void *params);
double M_delta_to_M_vir(double z, double mdelta, double delta);
double M_delta_to_M_vir_fast(double z, double mdelta, double delta);
double get_eq_for_M_delta(double x, void *params);
double M_vir_to_M_delta(double z, double mvir, double delta);

//integral with input tables
double tab_spline_and_integral(int Nbin, double *xlist, double *ylist, double *zlist);

//integral with input tables and include finite planes
double tab_spline_and_integral_forsim(int Nbin, double *xlist, double *ylist, double *zlist, int Nplane, double dplane, double *tablez);

//halo concentration by Diemer & Kravtosov 2015, https://arxiv.org/abs/1407.4730
double dlnP_dlnk(double lnk);
double nu_M(double z, double M);
double nu_M_fast(double z, double M);
double c_200c_DK15(double z, double M);
double c_vir_DK15_fast(double z, double M);
double get_eq_for_M_200c(double x, void *params);
double M_vir_to_M_200c(double z, double mvir);
