//input cosmological parameters
static double  Omega;
static double  OmegaLambda;
static double  w;
static double  Sigma8;
static double  HubbleParam;
const  double  OmegaBaryon =0.0456;
//unused parameters
const double  Gamma       =0.201;
const double  fnl         =0;
const double  delta_R     =2.441e-9;
const double  ns          =0.96;
const double  N_nu        =3.0;

const double  SEC_PER_YEAR=3.155e7;

const double T_CMB = 2.726; //K
const double freq_obs = 143.0; //GHz 

/* Cosmic abundance in plasma
element & fraction & mass fraction & fraction of free electrons
-----------------------------------------------------------------------------------------
H                 100           100                   100
He                8.5            34                    17
C, N, O, Ne     0.116          1.75                   0.9
Others          0.014          0.50                  0.23
-----------------------------------------------------------------------------------------
sum              108.63        136.25              118.13
*/
const double G0=4.299e-9; //km^2 Mpc/M_sun/s^2
const double mean_weight = 0.6; //136.25/(108.63+118.13);
const double Nx = 0.76; //primordial hydrogen abundance (P_e=(2+2*Nx)/(3+5*Nx)*P_gas)
const double Msun = 1.9884e30; //kg
const double kB = 6.94352394e-60; //km^2 Msun s^-2 K^-1
const double m_proton = 8.4118969e-58; //Msun

const double sigma_T = 0.665245854e-24; //(cm)^2 Tomson cross-section
const double m_elect = 0.510998902e6; //eV/c^2
const double planck_h = 4.135667e-15; //eV s

const double delta_h   = 200;
const double zhalo_max = 7.0;
const double Mhalo_max = 7.0e15;

//don't change below
const double  PI          =3.14159265358979323846; 
const double  GRAVITY     =6.672e-8;
const double  HUBBLE      =3.2407789e-18;   /* in h/sec */

const double  H0          =0.1;
const double  GCONST      =43007.1;
const double  Rlog_e10    =2.302585093; //ln10
const double  TINY        =1.0e-48;
const double  C           =2.998e5;
const double  R8=8000.0;    /* 8 Mpc/h */
const double  Mwidth      =0.05;
const int     NPOINTS= 20;
//                         Internally Defined Constants                       //
static const long double pi =  3.1415926535897932384626433832795029L; 
static const long double pi2 = 1.5707963267948966192313L;      // pi / 2
static const long double euler_gamma = 0.577215664901532860606512090L;
static const double auxiliary_asymptotic_cutoff = 48.0;
