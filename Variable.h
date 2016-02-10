/***************************************/
/***              EOS                ***/
/***       equation of state         ***/
/***          Stefan Typel           ***/
/***        new version 3.13         ***/
/***           2012/01/13            ***/
/***************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <map>
#include <vector>
/******************/
/***** macros *****/
/******************/

#define QUAD(a) ( (a) * (a) )
#define CUBE(a) ( (a) * (a) * (a) )
#define QUAR(a) ( (a) * (a) * (a) * (a) )

/*********************/
/***** constants *****/
/*********************/

#define N_AME03  3178  /* number of nuclei in AME03 + 2                  */
#define N_AME11  3288  /* number of nuclei in AME11                      */
#define N_SPIN   6800  /* number of nuclei in groundstate spin table     */

#define N_CL       11  /* number of clusters                             */
#define N_HYP       6  /* number of hyperons                             */
#define N_TMES      5  /* number of thermal mesons                       */
#define N_PART     31  /* N_CL+N_HYP+N_TMES+9 number of particles        */

#define N_MF        6  /* number of interacting meson fields             */

#define N_PH        4  /* maximum number of phases                       */

#define DIM_R     401  /* maximum number of meshpoints for radius (odd!) */

#define DIM_ACC 18446  /* 2*N_PART*DIM_R                                 */
#define REC_ACC     4  /* depth of recursion in acceleration + 1         */
#define DIM_LA      3  /* dimension of linear algebra calculations
= REC_ACC-1                                    */

#define DIM_N_B   302  /* dimension in n_b for radii table               */
#define DIM_Y_Q   101  /* dimension in y_q for radii table               */
#define A_NO_RANGE 350
#define WRITE_BLOCK_SIZE 100
#define WINDOW_SIZE 64

/***poke progress engine metronomically***/
#define POKE_P_E 1


/**********************/
/***** structures *****/
/**********************/
typedef struct
{
	int total_max_scal, max_no_of_procs;
	double total_exec_time, min_wall_time;
	int a_no[A_NO_RANGE];
	int max_scal[A_NO_RANGE];
	double min_time[A_NO_RANGE];
	int no_of_procs_to_be_allocated[A_NO_RANGE];
} data3;

typedef struct m2
{
	double a[WRITE_BLOCK_SIZE];
	double b[WRITE_BLOCK_SIZE];
	int c[WRITE_BLOCK_SIZE];
	int d[WRITE_BLOCK_SIZE];
	double e[WRITE_BLOCK_SIZE];
	double f[WRITE_BLOCK_SIZE];
	double g[WRITE_BLOCK_SIZE];
	double h[WRITE_BLOCK_SIZE];
	double i[WRITE_BLOCK_SIZE];
	double j[WRITE_BLOCK_SIZE];
	double k[WRITE_BLOCK_SIZE];
	double l[WRITE_BLOCK_SIZE];
} data2;

struct param
{
	int NMZ, T, NB;
};


struct field
{
	double m;  /* meson mass [1/fm]             */
	double m2; /* square of meson mass [1/fm^2] */
};

struct part
{
	int a;      /* baryon (mass) number                           */
	int n;      /* neutron number                                 */
	int z;      /* charge number                                  */
	int s;      /* strangeness number                             */
	int st;     /* statistic                                      */
	int in;     /* included in nse calculation                    */
	double g;   /* degeneracy factor                              */
	double tdgdt;/* T*temperature derivative of degeneracy factor */
	double g0;  /* original degeneracy factor for clusters        */
	double be;  /* binding energy [1/fm]                          */
	double m;   /* mass [1/fm]                                    */
	double fs;  /* coupling factor to sigma meson                 */
	double fo;  /* coupling factor to omega meson                 */
	double fr;  /* coupling factor to rho meson                   */
	double fd;  /* coupling factor to delta meson                 */
	double fp;  /* coupling factor to phi meson                   */
	double rms; /* rms radius of particle [fm]                    */
	double ag;  /* parameter in density distribution [1/fm]       */
	double mu;  /* chemical potential [1/fm]                      */
	double mm;  /* in medium mass [1/fm]                          */
};

struct solution
{
	int type;
	int conv;
	int nr;
	double mu_b;
	double mu_q;
	double mu_l;
	double f;
	double r;
};

struct kern
{
	int n;
	int z;
	int a;
	double b;
	double bea;
	double g0;
	double g;
	double m0;
	double m;
	double ms;
	double dmdo;
	double dmdr;
	double dmdt;
	double mu;
	double dlngdlnt;
};

/****************************/
/***** global variables *****/
/****************************/
MPI_Group new_group, orig_group;
MPI_Comm new_comm, old_comm;

MPI_Win win_q;
MPI_Win win_offs;

int lock=-2;

/* tmp vars for logging */
double gt000;
double wait_for_access;
/* -------------------- */

char newfile[1000];
int debug, cswitch;

static int count = 0;
int reset = 0;
MPI_File fh,fh1;

double arrr;
int disp_unit;
int pure_work_time;
char mytimefile[100], totaltimefile[100], foldername[100];

void resetcount()
{
	count = reset;
}

/* test */
double XP0, XP1, XP2, XP3, XXX;

/* control variables */
int NL, IREXT, IN_CL, IN_HYP, IN_TMES, XCL, METHOD, IN_R, IDX_TAB[3], DBDEP, XHYP;

/* numerical constants */
double PI, TPI, FPI, PI2, TPI2, RPI;

/* physical constants */
double HBARC, HBARC3, ALPHA, E2, G_G, AMU, M_PRO, M_NEU, CONV[3],
M_LAMBDA, M_SIGMAP, M_SIGMA0, M_SIGMAM, M_XI0, M_XIM,
M_PI0, M_PIP, M_K0, M_KP, M_ETA,
U_LAMBDA, U_SIGMA, U_XI, R_HYP, R_PHI;

/* parameters etc. */
double TT, N_B, Y_Q, Y_C, Y_S;

/* nonlinear RMF models */
double NL_GS, NL_GO, NL_GR, NL_G2, NL_G3, NL_C3, PARA_RMF;

/* density dependent RMF models */
double RHOREF, RHOSREF, COEFF[N_MF][6], LA_OMEGA, LA_RHO, CPL[N_MF][2], CPL2[N_MF][2];

/* interaction meson fields */
double MF[N_MF][DIM_R], MF_OLD[N_MF][DIM_R], S_MF[N_MF][DIM_R],
XMF[N_PH][N_MF][DIM_R];
struct field MESON[N_MF];

/* constituent particles */
int IN_PART[N_PART], IN_TON;
struct part NUCLEUS[N_AME03], PARTICLE[N_PART];

/* table of nuclei */
int N_NUC, NP_NUC;
double X_NUC[N_AME11], Y_NUC[N_AME11], O_TON, S_TON, G_TON, DM_TON, DDMDN_TON, DDMDT_TON;
struct kern NUC[N_AME11];

/* sorting */
int ISORT[N_AME11];
double VSORT[N_AME11];

/* self-energies */
double VV[N_PART][DIM_R], SS[N_PART][DIM_R], MM[N_PART][DIM_R],
VV_NUC[N_AME11], SS_NUC[N_AME11];

/* thermodynamical properties */
double E_DENS, F_DENS, G_DENS, O_DENS, S_DENS, PRES, F_ENT, F_ENE, F_PRE, C_DENS,
F_IE, F_IS, F_IG, F_IP, F_XE, F_XS, F_XG, F_XP;

/* phases and solutions */
int IN_PH[N_PH], IPHASE, IEXR[N_PH], ISOL;
double R_IN[DIM_N_B][DIM_Y_Q][N_PH];
struct solution SOL[N_PH];

/* densities */
double RHO_MF[N_MF][DIM_R], RHOC_MF[N_MF][DIM_R], RHO_PS[N_MF][DIM_R],
RHO_TOT[DIM_R], RHOI_TOT[DIM_R],
DENS[N_PART][DIM_R][3], DENSS[N_PART][DIM_R][3],
DENSC[N_PART][DIM_R][3], DENSCS[N_PART][DIM_R][3],
DENS_N[4][DIM_R], DENSS_N[4][DIM_R], DENS_FOLD[DIM_R],
DENS_TON[N_AME11], DENSS_TON[N_AME11], DENS_TON_B, DENS_TON_Q;

/* densities for saving */
double XDENS[N_PH][N_PART][DIM_R][3], XDENSS[N_PH][N_PART][DIM_R][3],
XDENSC[N_PH][N_PART][DIM_R][3], XDENSCS[N_PH][N_PART][DIM_R][3];

/* particle fractions */
double X_PART[N_PART], Y_PART[N_PART], X_H, Y_H;

/* nucleus */
double BEA_NUC, BFA_NUC, BSA_NUC, BCA_NUC;

/* WS cell */
int NR, NRP, NX;
double R_WS, V_WS, R1, R2, V1, V2, VX, AA_WS, ZZ_WS, NN_WS, A0_WS, Z0_WS, N0_WS,
AA_NUC, RMS[0],
DR, RVEC[3][DIM_R],
GS_S0[DIM_R], GS_SM[DIM_R], GS_SP[DIM_R],
GS_A0[N_MF][DIM_R], GS_AM[N_MF][DIM_R], GS_AP[N_MF][DIM_R],
GS_FP_P[DIM_R], GS_FP_0[DIM_R], GS_FP_M[DIM_R],
GS_FPP_P[DIM_R], GS_FPP_0[DIM_R], GS_FPP_M[DIM_R];

/* Coulomb correction */
double F_COUL1, F_COUL2, DMU;

/* chemical potentials */
double MU_BB, MU_QQ, MU_LL, F_MU_B, F_MU_Q, F_YQ, F_AA;

/* densities */
int F_ST, IBX;
double F_G, F_M, F_V, F_S, F_RHOP, F_RHOA, F_MD, F_MD2, F_KM, F_KP;

/* cluster shifts */
double SHIFT_F, SHIFT_FP, SHIFT_FT, DEG[N_CL], DEGT[N_CL], DEJ, DEJT,
DEPG[2][N_CL], DEPJ[3],
A_SC[N_CL], R_SC[N_CL], CL_A[N_CL], CL_B[N_CL],
CL_C[N_CL], CL_D[N_CL], CL_MU[N_CL], GGMM[N_CL], TDGGMM[N_CL],
VI[N_CL], TDVI[N_CL], BE_CL[N_CL], DBE_CL[4][N_CL], E_RES[N_CL];

/* convergence acceleration */
int M_ACC, N_ACC, M_ACC_MF, N_ACC_MF, M_ACC0, N_ACC0;
double VMB[DIM_ACC][REC_ACC], FMB[DIM_ACC][REC_ACC],
VMB0[DIM_ACC][REC_ACC], FMB0[DIM_ACC][REC_ACC],
DVMB[DIM_ACC][REC_ACC], DFMB[DIM_ACC][REC_ACC],
UMB[DIM_ACC][REC_ACC], GUMB[DIM_ACC][REC_ACC],
AMB[REC_ACC][REC_ACC], BMB[REC_ACC][REC_ACC], CMB[REC_ACC], GMB[REC_ACC],
VMF[DIM_R][REC_ACC], FMF[DIM_R][REC_ACC],
DVMF[DIM_R][REC_ACC], DFMF[DIM_R][REC_ACC],
UMF[DIM_R][REC_ACC], GUMF[DIM_R][REC_ACC],
AMF[DIM_R][REC_ACC], BMF[DIM_R][REC_ACC], CMF[REC_ACC], GMF[REC_ACC];

/* convergence checks */
int CI_NB, CI_YQ, CI_NBYQ, CI_RMF, CI_R, ICV_RMF;

/* matrix inversion */
double MAT[DIM_LA][DIM_LA], MAT2[DIM_LA][DIM_LA], MAT3[DIM_LA][DIM_LA],
XVEC[DIM_LA], YVEC[DIM_LA], YVEC2[DIM_LA], LUVV[DIM_LA];
int INDX[DIM_LA];
/*my file*/

FILE *myfile;

/* data files */
FILE *FILE_PARA, *FILE_IN, *FILE_OUT_THERMO, *FILE_OUT_COMPO, *FILE_OUT_EXTRA,
*FILE_R1, *FILE_R2, *FILE_R3, *FILE_PLOT, *FILE_PLOT2, *FILE_OUT_COMPOH;

/* temporary */
int IMFS[4], IN_C[12], ICTR[4], INDEX[REC_ACC], INDEX_MF[REC_ACC], ICOUNT;
