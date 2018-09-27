double PowerSpec(double k);
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

inline double efnn(double x, double);
inline double var4(double x, double);
inline double efn(double x, double);
inline double var3(double x, double);
inline double ddweight(double x, double);
inline double dweight(double x, double);
inline double weight(double x, double);
double dlogdsigma(double mass, double, double);

double sigma_m(double m, double *rsphere_return);
double fract_mass(double sig);
double sigdsigdr(double);
double dsig2dr2(double);
double unnsigma(double);
inline double evar2(double x, double);
inline double var2(double x, double);
double dndlogm(double logm);
double dndlogm_NG(double logm);
double dndlogm_MG(double logm);

void readCMB_and_do_spline();

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y);
void splie3(double x1a[], double x2a[], double x3a[], double ***ya, int l, int m, int n, double ***y2a);
void splin3(double x1a[], double x2a[], double x3a[], double ***ya, double ***y2a, int l, int m, int n, double x1, double x2, double x3, double *y);


double H_z(double z);
double H_z_1(double z,void *params);
double chi(double z);
double rho_nfw(double r,double z,double M);
double rho_s(double z,double M);
double r_s(double z,double M);
double r_vir(double z, double M);
double c_nfw(double z,double M);
double c_prada(double z, double M);
double delta_c(double z,double M);
double delta_v(double z);
double r_delta(double z, double Mass, double delta);
double get_eq_for_M_vir(double x, void *params);
double M_delta_to_M_vir(double z, double mdelta, double delta);

double get_eq_for_M_delta(double x, void *params);
double M_vir_to_M_delta(double z, double mvir, double delta);

double rho_s_fast(double z,double M);
double r_s_fast(double z,double M);
double c_nfw_fast(double z,double M);

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

double U_nfw(double M, void *params);
double U_nfw_fast(double log_ell, double logm, double z);

//2halo term
double halo_bias(double logm);
double halo_bias_shift(double logm);
double halo_bias_NG(double k, double logm, double z);
double halo_bias_MG(double logm);
double halo_bias_fast(double logm, double z);
double Window_halo_int(double logm, double z);
double Window_halo(double z);
double Window_halo_kdepend_int(double logm, double z);
double Window_halo_kdepend(double ell, double z);
double Window_halo_fast(double z);
double Window_halo_kdepend_fast(double ell, double z);
double Window_kappa(double z);

double Omega_de(double a);
double coeff1(double a);
double coeff2(double a);//u"+coeff1(a)u'+coeff2(a)u=0, u=D/a, '=d/dlna
double RungeKutta(double a_in,double a); //calculate linear density growth eq.

void stack_table_and_spline();

double chi_fast(double z);
void growth_spline();
double growth_for_any_w(double a);
double dndlogm_fast(double logm, double z);
double r_vir_fast(double logm, double z);

double Redshift(double x,double z_in);

//C_kappa
double P_kappa_int(double z);
double P_kappa(double ell);
double P_kappa_fast(double ell);
double P_kappa_discrete_lowres(double ell, double shell, double covd1, double covd2);
double P_kappa_discrete_highres(double ell, double shell, double covd1, double covd2, char *list);

//non-linear matter Pk
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

void stack_data_and_spline_Pk();

double dsig_dR_fast(double R);
double d2sig_dR2_fast(double R);
double sigma2_gauss_fast(double lnR);

//galaxy HOD
double HOD_gal(double M, int cen_or_sat);
double mean_ngal_int(double logM, double z);
double mean_ngal(double z);
double mean_ngal_fast(double z);

double dens_sat_int(double logM, double z);
double dens_sat(double z);

//galaxy auto power spectrum
double onehalo_gg_int(double logm, double k, double z, int opt);
double twohalo_gg_int(double logm, double k, double z, int opt);

double Pgg_1halo(double k, double z, int opt);
double Pgg_2halo(double k, double z, int opt);

//galaxy-matter auto power spectrum
double onehalo_gm_int(double logm, double k, double z, int opt);
double twohalo_gm_int(double logm, double k, double z, int opt);

double Pgm_1halo(double k, double z, int opt);
double Pgm_2halo(double k, double z, int opt);

//kappa power spectrum for comoving cordinate

void save_window_kappa(double z_l);
void save_window_kappa_per_k(double z_l, double calc_k);

double P_kappa_comv_int(double k, double z);
double P_kappa_comv(double k);

double P_kappa_comv_dSigma_int(double k, double z);
double P_kappa_comv_dSigma(double k);

double P_kappa_comv_gammat_int(double k, double z);
double P_kappa_comv_gammat(double k);

double window_kappa_dSigma(double z);
double window_kappa_gammat(double z);
double window_kappa_dSigma_fast(double z);
double window_kappa_gammat_fast(double z);

//void save_optimal_weight();
//double stack_lensing_optimal_weight_fast(double k, double zs);

//halo concentration by Diemer & Kravtosov 2015, https://arxiv.org/abs/1407.4730
double dlnP_dlnk(double lnk);
double nu_M(double z, double M);
double nu_M_fast(double z, double M);
double c_200c_DK15(double z, double M);
double c_vir_DK15_fast(double z, double M);
double get_eq_for_M_200c(double x, void *params);
double M_vir_to_M_200c(double z, double mvir);

void set_halo_conc_DK15(char *inputPk, double Omega_m, double Omega_b, double w0, double h0, double ns_index);
void free_halo_conc_DK15();
