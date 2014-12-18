/** Old Version!!! Needs to be checked **/

/**************************************/
/***              EOS               ***/
/***       equation of state        ***/
/***          Stefan Typel          ***/
/***        new version 3.13        ***/
/***           2012/01/13           ***/
/**************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
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

/**********************/
/***** structures *****/
/**********************/

typedef struct m2
{
    double a[175];
    double b[175];
    int c[175];
    int d[175];
    double e[175];
    double f[175];
    double g[175];
    double h[175];
    double i[175];
    double j[175];
    double k[175];
    double l[175];

} data2;

struct param
{
    int NMZ , T , NB;
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

int debug;

static int count=0;
int reset=0;

void resetcount()
{
    count=reset;
}


/* test */
double XP0,XP1,XP2,XP3,XXX;

/* control variables */
int NL,IREXT,IN_CL,IN_HYP,IN_TMES,XCL,METHOD,IN_R,IDX_TAB[3],DBDEP,XHYP;

/* numerical constants */
double PI,TPI,FPI,PI2,TPI2,RPI;

/* physical constants */
double HBARC,HBARC3,ALPHA,E2,G_G,AMU,M_PRO,M_NEU,CONV[3],
       M_LAMBDA,M_SIGMAP,M_SIGMA0,M_SIGMAM,M_XI0,M_XIM,
       M_PI0,M_PIP,M_K0,M_KP,M_ETA,
       U_LAMBDA,U_SIGMA,U_XI,R_HYP,R_PHI;

/* parameters etc. */
double TT,N_B,Y_Q,Y_C,Y_S;

/* nonlinear RMF models */
double NL_GS,NL_GO,NL_GR,NL_G2,NL_G3,NL_C3,PARA_RMF;

/* density dependent RMF models */
double RHOREF,RHOSREF,COEFF[N_MF][6],LA_OMEGA,LA_RHO,CPL[N_MF][2],CPL2[N_MF][2];

/* interaction meson fields */
double MF[N_MF][DIM_R],MF_OLD[N_MF][DIM_R],S_MF[N_MF][DIM_R],
       XMF[N_PH][N_MF][DIM_R];
struct field MESON[N_MF];

/* constituent particles */
int IN_PART[N_PART],IN_TON;
struct part NUCLEUS[N_AME03],PARTICLE[N_PART];

/* table of nuclei */
int N_NUC,NP_NUC;
double X_NUC[N_AME11],Y_NUC[N_AME11],O_TON,S_TON,G_TON,DM_TON,DDMDN_TON,DDMDT_TON;
struct kern NUC[N_AME11];

/* sorting */
int ISORT[N_AME11];
double VSORT[N_AME11];

/* self-energies */
double VV[N_PART][DIM_R],SS[N_PART][DIM_R],MM[N_PART][DIM_R],
       VV_NUC[N_AME11],SS_NUC[N_AME11];

/* thermodynamical properties */
double E_DENS,F_DENS,G_DENS,O_DENS,S_DENS,PRES,F_ENT,F_ENE,F_PRE,C_DENS,
       F_IE,F_IS,F_IG,F_IP,F_XE,F_XS,F_XG,F_XP;

/* phases and solutions */
int IN_PH[N_PH],IPHASE,IEXR[N_PH],ISOL;
double R_IN[DIM_N_B][DIM_Y_Q][N_PH];
struct solution SOL[N_PH];

/* densities */
double RHO_MF[N_MF][DIM_R],RHOC_MF[N_MF][DIM_R],RHO_PS[N_MF][DIM_R],
       RHO_TOT[DIM_R],RHOI_TOT[DIM_R],
       DENS[N_PART][DIM_R][3],DENSS[N_PART][DIM_R][3],
       DENSC[N_PART][DIM_R][3],DENSCS[N_PART][DIM_R][3],
       DENS_N[4][DIM_R],DENSS_N[4][DIM_R],DENS_FOLD[DIM_R],
       DENS_TON[N_AME11],DENSS_TON[N_AME11],DENS_TON_B,DENS_TON_Q;

/* densities for saving */
double XDENS[N_PH][N_PART][DIM_R][3],XDENSS[N_PH][N_PART][DIM_R][3],
       XDENSC[N_PH][N_PART][DIM_R][3],XDENSCS[N_PH][N_PART][DIM_R][3];

/* particle fractions */
double X_PART[N_PART],Y_PART[N_PART],X_H,Y_H;

/* nucleus */
double BEA_NUC,BFA_NUC,BSA_NUC,BCA_NUC;

/* WS cell */
int NR,NRP,NX;
double R_WS,V_WS,R1,R2,V1,V2,VX,AA_WS,ZZ_WS,NN_WS,A0_WS,Z0_WS,N0_WS,
       AA_NUC,RMS[0],
       DR,RVEC[3][DIM_R],
       GS_S0[DIM_R],GS_SM[DIM_R],GS_SP[DIM_R],
       GS_A0[N_MF][DIM_R],GS_AM[N_MF][DIM_R],GS_AP[N_MF][DIM_R],
       GS_FP_P[DIM_R],GS_FP_0[DIM_R],GS_FP_M[DIM_R],
       GS_FPP_P[DIM_R],GS_FPP_0[DIM_R],GS_FPP_M[DIM_R];

/* Coulomb correction */
double F_COUL1,F_COUL2,DMU;

/* chemical potentials */
double MU_BB,MU_QQ,MU_LL,F_MU_B,F_MU_Q,F_YQ,F_AA;

/* densities */
int F_ST,IBX;
double F_G,F_M,F_V,F_S,F_RHOP,F_RHOA,F_MD,F_MD2,F_KM,F_KP;

/* cluster shifts */
double SHIFT_F,SHIFT_FP,SHIFT_FT,DEG[N_CL],DEGT[N_CL],DEJ,DEJT,
       DEPG[2][N_CL],DEPJ[3],
       A_SC[N_CL],R_SC[N_CL],CL_A[N_CL],CL_B[N_CL],
       CL_C[N_CL],CL_D[N_CL],CL_MU[N_CL],GGMM[N_CL],TDGGMM[N_CL],
       VI[N_CL],TDVI[N_CL],BE_CL[N_CL],DBE_CL[4][N_CL],E_RES[N_CL];

/* convergence acceleration */
int M_ACC,N_ACC,M_ACC_MF,N_ACC_MF,M_ACC0,N_ACC0;
double VMB[DIM_ACC][REC_ACC],FMB[DIM_ACC][REC_ACC],
       VMB0[DIM_ACC][REC_ACC],FMB0[DIM_ACC][REC_ACC],
       DVMB[DIM_ACC][REC_ACC],DFMB[DIM_ACC][REC_ACC],
       UMB[DIM_ACC][REC_ACC],GUMB[DIM_ACC][REC_ACC],
       AMB[REC_ACC][REC_ACC],BMB[REC_ACC][REC_ACC],CMB[REC_ACC],GMB[REC_ACC],
       VMF[DIM_R][REC_ACC],FMF[DIM_R][REC_ACC],
       DVMF[DIM_R][REC_ACC],DFMF[DIM_R][REC_ACC],
       UMF[DIM_R][REC_ACC],GUMF[DIM_R][REC_ACC],
       AMF[DIM_R][REC_ACC],BMF[DIM_R][REC_ACC],CMF[REC_ACC],GMF[REC_ACC];

/* convergence checks */
int CI_NB,CI_YQ,CI_NBYQ,CI_RMF,CI_R,ICV_RMF;

/* matrix inversion */
double MAT[DIM_LA][DIM_LA],MAT2[DIM_LA][DIM_LA],MAT3[DIM_LA][DIM_LA],
       XVEC[DIM_LA],YVEC[DIM_LA],YVEC2[DIM_LA],LUVV[DIM_LA];
int INDX[DIM_LA];
/*my file*/

FILE *myfile;

/* data files */
FILE *FILE_PARA,*FILE_IN,*FILE_OUT_THERMO,*FILE_OUT_COMPO,*FILE_OUT_EXTRA,
     *FILE_R1,*FILE_R2,*FILE_R3,*FILE_PLOT,*FILE_PLOT2,*FILE_OUT_COMPOH;

/* temporary */
int IMFS[4],IN_C[12],ICTR[4],INDEX[REC_ACC],INDEX_MF[REC_ACC],ICOUNT;

/*******************************/
/***** function prototypes *****/
/*******************************/

/* basics routines */
void worker (int mpitask_id, long int MAX_SIZE, int task[]);
void manager(int number_of_processors, int argc , char * argv[]);
void writer(void);

data2* myfunction(char *,int,int,int,int,data2*);
char * changestring(char * ,char*,char *);
void call_error(int);

/* initialization */
int init_eos(int,int,int,char*),init_rmf(int),init_ton(int,char*),init_nuclei_ame03(int,char*),
    init_parameters(void),init_particles(int),init_dens(void),
    init_self_energies(void),init_mf(void);

/* coupling functions */
int get_cpl(double,int);

/* nuclei table */
data2* get_nuc_table(int,int,int,int,int,int,int *,int,int,int,int,int,char *,data2 *,int);

int  get_nuc_table_old(int,int,int,int,int,int,int *,int,int,char *),
     fit_radius_nuc(int,double),fit_nmz_nuc(int,double),
     get_nucleus(double,double,char*),get_ton_mass(double,double),get_prop_ton(void),
     get_ton_mass_shift(int,double,double),get_ton_g(int);
double get_nucleus_old(double,double);

/* EoS table */
int get_eos_table(int,int,int,int,int,int,int *,int,int,int,int *,int,int,char *),
    get_parameter_tables(int,int,int,double,double,double,double,char *);

/* basic EoS routines */
int solve_rmf(int,double),cycle_solve(int,double,double),
    cycle_solve_old(int,double,double),
    get_eos_composition(int,char*),
    get_properties(void),
    get_eos_local(double,double,double,int,int,int,int,int,int,
                  int *,int,int,int,int,int,double,double,char*);
double get_eos_local_beta(double,double,int,int,int,int,int,int,
                          int *,int,int,int,int,int,double,double,char*);

/* meson fields and self-energies */
int get_mf(int),get_self_energies(int),get_sources(void),
    get_mf_greens(int,double *);

/* chemical potentials*/
double get_chemical_potential_lep(void);
int get_chemical_potential_bar(double,double),
    fit_asymmetry(double,double,double,double),
    fit_density(double,double,double);

/* properties */
int get_prop_part(int),get_prop(double),get_anz_nuc(void);

/* densities */
double get_rho_part(double,int,int),get_n_lep(double,double,int),
       disp_rel(double,int),inv_disp_rel(double,int),
       get_simint(int,double,int,double,double,double,int);
int get_rho(double),get_rhos(double),get_simint_bounds_f(double,int),
    get_simint_bounds_b(double,int),get_nb_yq(double,double,int),
    get_dens_ton(int);

/* cluster shifts */
int get_shift(double,double,double,double),get_dbe_p_g(int),get_dbe_p_j(void),
    init_be_cl(void),get_be_cl3(int,double,double);

/* WS cell */
int fit_radius2(int,double),fit_aa4_new(int,double),init_discr_ws_cell(int),
    get_results_ws_cell(char *),fit_aa4(int,double),fit_aa4_old(int,double),
    get_coul_corr(void);
double solve_a(int,int,double);

/* solutions */
int init_sol(void),iexs_sol(int),save_sol(int),select_sol(int),
    get_fractions(void);

/* convergence acceleration */
int init_acc(void),get_acc(void),init_acc_mf(int),get_acc_mf(int),
    init_acc0(void),get_acc0(void);

/* matrix inversion */
int get_mat_inv(int),lu_decomposition(int),lu_backsubstitution(int),
    improve(int);

/* functions */
int fold_dens(int);
double f_ws(double,double,double),r_integrate(int,int,double *),
       trapez(double *),func_fold(double,double,double),
       simpson_mod(int,double *),f_func(double,int,double,int),
       func_fermi(double,double,int),func_bose(double,double,int),
       get_k0_red(double),get_k1_red(double),get_k2_red(double),experfc(double);

/* sorting */
int heapsort(int);

/* test */
int get_test(void),get_test2(void);

/*********************/
/***** functions *****/
/*********************/

/*****************************************************************************/
//data2* myfunction(char *s,int t_nmz , int t_nb , int t_t, int a_0, data2*point)
//{
//    int ipara,inuc,imb,inse,in_g,in_e,in_mu,in_tau,in_nu,
//        in_cl[N_CL],in_hyp,in_tmes,in_ton,in_phase[N_PH],in_r,in_tab;
//
//    XXX = 0.;
//
//    /***** INPUT *****/
//
//    /* RMF parametrization */
//    /* 0: no interaction, 2: DD2, 3: TM1, 4: DD-MEdelta, default: DD */
//    ipara = 2;
//
//    /* calculation of nuclei or EoS */
//    /* 1: nuclei */
//    /* else: EoS*/
//    inuc = 1;
//
//    /* for calculation of EoS */
//    /* 1: NSE calculation, default: RMF */
//    /* noch gebraucht? in_ton reicht */
//    inse = 0;
//
//    /* for light clusters */
//    /* 1: use non-relativistic Maxwell-Boltzmann statistics */
//    /* 2: use relativistic Maxwell-Boltzmann statistics */
//    /* else: use original statistics */
//    imb = 1;
//
//    /* include photon */
//    in_g = 0;
//    /* include electron */
//    in_e = 0;
//    /* include muon */
//    in_mu = 0;
//    /* include tauon */
//    in_tau = 0;
//    /* include neutrinos */
//    in_nu = 0;
//    /* include deuteron */
//    in_cl[0] = 0;
//    /* include triton */
//    in_cl[1] = 0;
//    /* include helion */
//    in_cl[2] = 0;
//    /* include alpha */
//    in_cl[3] = 0;
//    /* include np 3S1 channel */
//    in_cl[4] = 0;
//    /* include np 1S0 channel */
//    in_cl[5] = 0;
//    /* include nn 1S0 channel */
//    in_cl[6] = 0;
//    /* include pp 1S0 channel */
//    in_cl[7] = 0;
//    /* include triton continuum */
//    in_cl[8] = 0;
//    /* include helion continuum */
//    in_cl[9] = 0;
//    /* include alpha continuum */
//    in_cl[10] = 0;
//    /* include hyperons */
//    in_hyp = 0;
//    /* include thermal mesons */
//    in_tmes = 0;
//    /* include table of nuclei (only for eos calculation) */
//    in_ton = 0;
//    if (inuc == 1) in_ton = 0;
//
//    /* include drop phase */
//    in_phase[1] = 0;
//    /* include bubble phase */
//    in_phase[2] = 0;
//    /* include hole phase */
//    in_phase[3] = 0;
//
//    /* 1: read data2files for radii */
//    in_r = 0;
//
//    /* 1: standard input file eos_table.in
//       else: input file eos_tab.in */
//    in_tab = 1;
//
//    /* 0: dependence of energy shifts on meson fields
//       else: dependence of energy shifts on densities */
//    DBDEP = 0;
//    int kk=0;
//    /***** INITIALIZATION *****/
//
//    init_eos(ipara,imb,in_ton,s);
//
//    /***** CALCULATION *****/
//
//    if (inuc == 1)
//    {
//        point=get_nuc_table(inse,in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,t_nmz,t_nb,t_t,s,point, a_0);
//    }
//    else
//    {
//        get_eos_table(inse,in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
//                      in_ton,in_phase,in_r,in_tab,s);
//    }
//
//    /*    for(kk=0; kk<count; kk++)
//            fprintf(myfile,"%e %e %i %i %e %e %e %e %e %e %e %e\n",point->a[kk],point->b[kk],point->c[kk], point->d[kk],point->e[kk],point->f[kk],point->g[kk],point->h[kk],point->i[kk],point->j[kk],point->k[kk], point->l[kk]) ;
//    */
//    return point;
//
//}
/*****************************************************************************/
data2* get_nuc_table(int inse,int in_g,int in_e,int in_mu,int in_tau,int in_nu,
                     int *in_cl,int in_hyp,int in_tmes,int t_nmz, int t_nb, int t_t, char *cnsfile,data2* point, int a_0 )
{
    int iwr,itmp,icl,iph,in_cell,inuc,in_a,nmz_min,nmz_max,aa_max,nmz,ia,
        ic,ic2,in_ton,iz,shift;
    /*float tmp;*/
    double t,n_b,y_q,a,z,dnmz,daa_max,tmp1,tmp2;

    iwr = 0;

    IN_R = 0;
    inse = 0;
    in_g = in_mu = in_tau = in_nu = in_hyp = in_tmes = in_ton = 0;
    for (icl=0; icl<N_CL; icl++) in_cl[icl] = 0;
    IN_PH[0] = IN_PH[1] = 1;
    for (iph=2; iph<N_PH; iph++) IN_PH[iph] = 0;
    in_cell = 1;
    inuc = 1;

    /* read data file nuc.in */
    //FILE_IN = fopen(s,"r");
    //fscanf(FILE_IN,"%i", &itmp);
    nmz_min = t_nmz;

    //fscanf(FILE_IN,"%i", &itmp);
    nmz_max = t_nmz;

    /*nmz_max = nmz_min;*/
    /*fscanf(FILE_IN,"%i", &itmp);
      aa_max = itmp;*/
    aa_max = 350;
    daa_max = (double)aa_max+0.5;
    if (1 == 0)
    {
        if(debug==1)if(debug==1)fprintf(myfile," T and n_b not read from table\n");
        t = 0.;
        n_b = 0.000;
        if(debug==1)fprintf(myfile," T = %e MeV\n", t);
        if(debug==1)fprintf(myfile," n_b = %e fm^-3\n", n_b);
        shift = (int)(n_b*4000.+0.001);
        if(debug==1)fprintf(myfile," shift = %i\n", shift);
    }
    else
    {
        /*
        fscanf(FILE_IN,"%f", &tmp);
        t = (double)tmp;
        fscanf(FILE_IN,"%f", &tmp);
        n_b = (double)tmp;
        */
        //fscanf(FILE_IN,"%i", &itmp);
        t = (double)t_t/10.;
        //fscanf(FILE_IN,"%i", &itmp);
        n_b = (double)t_nb/1000.;
        shift = t_nb*4;
        if(debug==1)fprintf(myfile," shift = %i\n", shift);
    }
    //fclose(FILE_IN);

    //FILE_OUT_EXTRA  = fopen(changestring("nuc",cnsfile,".dat"),"w");
    if(debug==1)
        FILE_PLOT       = fopen(changestring("plot",cnsfile,".dat"),"w");
    if(debug==1)FILE_PLOT2      = fopen(changestring("plot2",cnsfile,".dat"),"w");


    /* non-existence of radius solutions before first run */
    for (iph=0; iph<N_PH; iph++)
    {
        IEXR[iph] = 0;
        IMFS[iph] = 0;
    }

    /*
    if (n_a > 1) {
      da = (a_max-a_min)/(double)(n_a-1);
    }
    else {
      da = 0.;
    }
    */

    /*
    if (nmz_min == -20) {
      fprintf(FILE_OUT_EXTRA,"# %e\n", n_b);
      fprintf(FILE_OUT_EXTRA,"# %e\n", t);
    }
    */


    if (1 == 0)
    {
        if (1 == 1)
        {
            /* only single nucleus */

            ia = 60;
            iz = 30;

            n_b = 0.;
            t = 0.;

            if(debug==1)fprintf(myfile," T = %e MeV\n", t);
            if(debug==1)fprintf(myfile," n = %e fm^-3\n", n_b);

            if (n_b > 0.)
            {
                in_e = 1;
            }
            else
            {
                in_e = 0;
            }

            nmz = ia-2*iz;
            a = (double)ia;
            z = (double)iz;
            y_q = 0.5*(1.-(double)nmz/a);
            /***************************************************************/
            get_eos_local(t,n_b,y_q,inse,
                          in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                          in_cell,in_ton,inuc,a,z,cnsfile);
            /***************************************************************/
            fprintf(FILE_OUT_EXTRA," %i %i %e %e %e %e %e %e %e %e %e %e\n",
                    nmz, ia, t, N0_WS, Z0_WS, R_WS,
                    BEA_NUC*HBARC, BFA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                    RMS[0], RMS[1]);
            printf(" %f %f\n", BEA_NUC*HBARC, BEA_NUC*HBARC*(N0_WS+Z0_WS));
            get_results_ws_cell(cnsfile);
        }
        else
        {

            /* find nucleus with maximum binding energy */
            n_b = 0.;
            in_e = 0;
            /* starting point */
            ia = 55;
            iz = 22;

            nmz = ia-2*iz;
            a = (double)ia;
            z = (double)iz;
            y_q = 0.5*(1.-(double)nmz/a);
            /***************************************************************/
            get_eos_local(t,n_b,y_q,inse,
                          in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                          in_cell,in_ton,inuc,a,z,cnsfile);
            /***************************************************************/
            /*printf(" A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);*/
            do
            {
                tmp1 = BEA_NUC;
                ia += 1;
                iz +=1;
                nmz = ia-2*iz;
                a = (double)ia;
                z = (double)iz;
                y_q = 0.5*(1.-(double)nmz/a);
                /***************************************************************/
                get_eos_local(t,n_b,y_q,inse,
                              in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                              in_cell,in_ton,inuc,a,z,cnsfile);
                /***************************************************************/
                /*printf(" A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);*/
            }
            while (BEA_NUC > tmp1);
            ia -= 1;
            iz -= 1;
            BEA_NUC = tmp1;
            do
            {
                tmp1 = BEA_NUC;
                ia -= 1;
                iz -=1;
                nmz = ia-2*iz;
                a = (double)ia;
                z = (double)iz;
                y_q = 0.5*(1.-(double)nmz/a);
                /***************************************************************/
                get_eos_local(t,n_b,y_q,inse,
                              in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                              in_cell,in_ton,inuc,a,z,cnsfile);
                /***************************************************************/
                /*printf(" A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);*/
            }
            while (BEA_NUC > tmp1);
            ia += 1;
            iz += 1;
            BEA_NUC = tmp1;

            do
            {
                if(debug==1)fprintf(myfile," A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);
                tmp2 = BEA_NUC;
                ia += 1;
                nmz = ia-2*iz;
                a = (double)ia;
                z = (double)iz;
                y_q = 0.5*(1.-(double)nmz/a);
                /***************************************************************/
                get_eos_local(t,n_b,y_q,inse,
                              in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                              in_cell,in_ton,inuc,a,z,cnsfile);
                /***************************************************************/
                /*printf(" A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);*/
                do
                {
                    tmp1 = BEA_NUC;
                    ia += 1;
                    iz +=1;
                    nmz = ia-2*iz;
                    a = (double)ia;
                    z = (double)iz;
                    y_q = 0.5*(1.-(double)nmz/a);
                    /***************************************************************/
                    get_eos_local(t,n_b,y_q,inse,
                                  in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                                  in_cell,in_ton,inuc,a,z,cnsfile);
                    /***************************************************************/
                    /*printf(" A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);*/
                }
                while (BEA_NUC > tmp1);
                ia -= 1;
                iz -= 1;
                BEA_NUC = tmp1;
                do
                {
                    tmp1 = BEA_NUC;
                    ia -= 1;
                    iz -=1;
                    nmz = ia-2*iz;
                    a = (double)ia;
                    z = (double)iz;
                    y_q = 0.5*(1.-(double)nmz/a);
                    /***************************************************************/
                    get_eos_local(t,n_b,y_q,inse,
                                  in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                                  in_cell,in_ton,inuc,a,z,cnsfile);
                    /***************************************************************/
                    /*printf(" A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);*/
                }
                while (BEA_NUC > tmp1);
                ia += 1;
                iz += 1;
                BEA_NUC = tmp1;
            }
            while (BEA_NUC > tmp2);
            ia -= 1;
            BEA_NUC = tmp2;


            if(debug==1)fprintf(myfile," A Z  BE/A %i %i  %e\n", ia, iz, BEA_NUC*HBARC);





        }
    }
    else
    {

        if(debug==1)fprintf(myfile," T = %e MeV\n", t);
        if(debug==1)fprintf(myfile," n = %e fm^-3\n", n_b);

        if (1 == 0)   /* selected nuclei ??? */
        {
            /*printf(" %i %i\n", nmz_min, nmz_max);*/
            ic2 = 0;
            nmz = nmz_min;
            if(debug==1)fprintf(myfile," n-z = %i\n", nmz);
            ic = 1;
            for (in_a=0; in_a<3; in_a++)
            {
                ia = nmz_max+4-2*in_a;
                if(debug==1)fprintf(myfile," a = %i\n", ia);
                if ((nmz+ia)%2 == 1)
                {
                    if(debug==1)fprintf(myfile," wrong input\n");
                    //exit(0);
                }
                a = (double)ia;
                y_q = 0.5*(1.-(double)nmz/a);
                z = a*y_q;
                if(debug==1)fprintf(myfile," a z ic %f %f %i\n", a, z, ic);
                if (n_b > 0.)
                {
                    in_e = 1;
                }
                else
                {
                    in_e = 0;
                }
                if (ic == 1)
                {
                    /***************************************************************/
                    if((int)get_eos_local(t,n_b,y_q,inse,
                                          in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                                          in_cell,in_ton,inuc,a,z,cnsfile)==999)
                        continue;
                    /***************************************************************/
                    if (n_b == 0.)
                    {
                        /*printf(" %e %e %e %i %i\n", A0_WS, Z0_WS, N0_WS, ISOL, ICV_RMF);*/
                        /* no plot beyond dripline for zero density calculation */
                        if (((a-A0_WS) < 1.e-06) && (Z0_WS > 2.) && (N0_WS > 2.)
                                && (ISOL == 1) && (ICV_RMF == 1))
                        {
                            if (A0_WS < daa_max)
                            {
                                fprintf(FILE_OUT_EXTRA," %e %e %i %i %e %e %e %e %e %e %e %e\n",
                                        t, n_b, nmz, ia, N0_WS, Z0_WS,
                                        BEA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                                        RMS[0], RMS[1], R_WS);
                                /*fprintf(FILE_OUT_EXTRA," %i %i %e %e %e %e %e %e %e %e %e\n",
                                  nmz, ia, n_b, N0_WS, Z0_WS, R_WS,
                                  BEA_NUC*HBARC, BFA_NUC*HBARC, BSA_NUC,
                                  RMS[0], RMS[1]);*/
                                ic2 += 1;
                            }
                        }
                        else
                        {
                            ic = 0;
                        }
                    }
                    else
                    {
                        if ((Z0_WS > 2.) && (N0_WS > 2.))
                        {
                            if (A0_WS < daa_max)
                            {
                                fprintf(FILE_OUT_EXTRA," %e %e %i %i %e %e %e %e %e %e %e %e\n",
                                        t, n_b, nmz, ia, N0_WS, Z0_WS,
                                        BEA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                                        RMS[0], RMS[1], R_WS);
                                ic2 += 1;
                            }
                        }
                        else
                        {
                            ic = 0;
                        }
                    }
                }
            }
        }
        else
        {

            ic2 = 0;
            for (nmz=nmz_min; nmz<(nmz_max+1); nmz++)
            {
                if(debug==1)fprintf(myfile," n-z = %i\n", nmz);
                dnmz = (double)nmz;
                /* starting point */
                /* a_0 = 80+shift;
                 if (nmz > 30) a_0 += 2*nmz;
                 a_0 = 80;
                */              // if (a_0 > aa_max) a_0 = aa_max;
                //if ((nmz+100)%2 == 1) a_0 += 1;

                ic = 1;
                if ( (nmz<a_0) && (nmz> (-1)*a_0))
                {
                    /*for (in_a=0; in_a<(a_0/2); in_a++)
                    {*/
                    ia = a_0;//-2*in_a;
                    a = (double)ia;
                    y_q = 0.5*(1.-(double)nmz/a);
                    z = a*y_q;
                    /*printf(" a z ic %f %f %i\n", a, z, ic);*/
                    if (n_b > 0.)
                    {
                        in_e = 1;
                    }
                    else
                    {
                        in_e = 0;
                    }
                    if (ic == 1)
                    {

                        /***************************************************************/
                        if((int)get_eos_local(t,n_b,y_q,inse,in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,in_cell,in_ton,inuc,a,z,cnsfile)==999)
                            return point;//continue; ///SHOULD BE "return point;" HERE??
                        /***************************************************************/

                        if (n_b == 0.)
                        {
                            /*printf(" %e %e %e %i %i\n", A0_WS, Z0_WS, N0_WS, ISOL, ICV_RMF);*/
                            /* no plot beyond dripline for zero density calculation */
                            if (((a-A0_WS) < 1.e-06) && (Z0_WS > 2.) && (N0_WS > 2.)
                                    && (ISOL == 1) && (ICV_RMF == 1))
                            {
                                if (A0_WS < daa_max)
                                {
                                    //	  fprintf(FILE_OUT_EXTRA," %e %e %i %i %e %e %e %e %e %e %e %e\n",
                                    //		  t, n_b, nmz, ia, N0_WS, Z0_WS,
                                    //		  BEA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                                    //		  RMS[0], RMS[1], R_WS);

                                    point->a[count]=t;
                                    point->b[count]=n_b;
                                    point->c[count]=nmz;
                                    point->d[count]=ia;
                                    point->e[count]=N0_WS;
                                    point->f[count]=Z0_WS;
                                    point->g[count]=BEA_NUC*HBARC;
                                    point->h[count]=BSA_NUC;
                                    point->i[count]=BCA_NUC*HBARC;
                                    point->j[count]=RMS[0];
                                    point->k[count]=RMS[1];
                                    point->l[count]=R_WS;
                                    count++;

                                    //		fprintf(myfile,"%e %e %i %i %e %e %e %e %e %e %e %e\n",point->a[count],point->b[count],point->c[count]
                                    //		,point->d[count],point->e[count],point->f[count],point->g[count],point->h[count],point->i[count],point->j[count],point->k[count],
                                    //		point->l[count])  ;



                                    /*fprintf(FILE_OUT_EXTRA," %i %i %e %e %e %e %e %e %e %e %e\n",
                                      nmz, ia, n_b, N0_WS, Z0_WS, R_WS,
                                      BEA_NUC*HBARC, BFA_NUC*HBARC, BSA_NUC,
                                      RMS[0], RMS[1]);*/
                                    ic2 += 1;
                                }
                            }
                            else
                            {
                                ic = 0;
                            }
                        }
                        else
                        {
                            if ((Z0_WS > 2.) && (N0_WS > 2.))
                            {
                                if (A0_WS < daa_max)
                                {
                                    //		fprintf(FILE_OUT_EXTRA," %e %e %i %i %e %e %e %e %e %e %e %e\n",
                                    //		  t, n_b, nmz, ia, N0_WS, Z0_WS,
                                    //		 BEA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                                    //		  RMS[0], RMS[1], R_WS);

                                    point->a[count]=t;
                                    point->b[count]=n_b;
                                    point->c[count]=nmz;
                                    point->d[count]=ia;
                                    point->e[count]=N0_WS;
                                    point->f[count]=Z0_WS;
                                    point->g[count]=BEA_NUC*HBARC;
                                    point->h[count]=BSA_NUC;
                                    point->i[count]=BCA_NUC*HBARC;
                                    point->j[count]=RMS[0];
                                    point->k[count]=RMS[1];
                                    point->l[count]=R_WS;
                                    count++;


                                    //		fprintf(myfile,"%e %e %i %i %e %e %e %e %e %e %e %e\n",point->a[count],point->b[count],point->c[count]
                                    //		,point->d[count],point->e[count],point->f[count],point->g[count],point->h[count],point->i[count],point->j[count],point->k[count],
                                    //		point->l[count])  ;


                                    ic2 += 1;
                                }
                            }
                            else
                            {
                                ic = 0;
                            }
                        }
                    }
                }
            }
            ic = 1;
            /*for (in_a=1; in_a<((aa_max-a_0)/2+1); in_a++)
            {
                ia = a_0+2*in_a;
                a = (double)ia;
                y_q = 0.5*(1.-(double)nmz/a);
                z = a*y_q;
                if(debug==1)fprintf(myfile," a z ic %f %f %i\n", a, z, ic);
                if (n_b > 0.)
                {
                    in_e = 1;
                }
                else
                {
                    in_e = 0;
                }
                if (ic == 1) // 19/03/2014
                {
                    itime = MPI_Wtime();
                    /***************************************************************/
            /*   if((int)get_eos_local(t,n_b,y_q,inse,
                                     in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                                     in_cell,in_ton,inuc,a,z,cnsfile)==999)
                   continue;
               /***************************************************************/


            /*
                                    ftime = MPI_Wtime()-itime;
                                if (n_b == 0.)
                                {
                                    /* no plot beyond dripline for zero density calculation */
            /*   if (((a-A0_WS) < 1.e-06) && (Z0_WS > 2.) && (N0_WS > 2.) && (ICV_RMF == 1)
                       && (ia <= aa_max) && (ISOL == 1))
               {
                   if (A0_WS < daa_max)
                   {
                       //	 fprintf(FILE_OUT_EXTRA," %e %e %i %i %e %e %e %e %e %e %e %e\n",
                       //		  t, n_b, nmz, ia, N0_WS, Z0_WS,
                       //		  BEA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                       //		  RMS[0], RMS[1], R_WS);
                       point->a[count]=t;
                       point->b[count]=n_b;
                       point->c[count]=nmz;
                       point->d[count]=ia;
                       point->e[count]=N0_WS;
                       point->f[count]=Z0_WS;
                       point->g[count]=BEA_NUC*HBARC;
                       point->h[count]=BSA_NUC;
                       point->i[count]=BCA_NUC*HBARC;
                       point->j[count]=RMS[0];
                       point->k[count]=RMS[1];
                       point->l[count]=R_WS;
                           point->m[count]=ftime;


                       //		fprintf(myfile,"%e %e %i %i %e %e %e %e %e %e %e %e\n",point->a[count],point->b[count],point->c[count]
                       //		,point->d[count],point->e[count],point->f[count],point->g[count],point->h[count],point->i[count],point->j[count],point->k[count],
                       //		point->l[count])  ;


                       count++;







                       ic2 =+ 1;
                   }
               }
               else
               {
                   ic = 0;
               }
            }
            else
            {
               if ((Z0_WS > 2.) && (N0_WS > 2.) && (ia <= aa_max))
               {
                   if (A0_WS < daa_max)
                   {
                       //	   fprintf(FILE_OUT_EXTRA," %e %e %i %i %e %e %e %e %e %e %e %e\n",
                       //		  t, n_b, nmz, ia, N0_WS, Z0_WS,
                       //		  BEA_NUC*HBARC, BSA_NUC, BCA_NUC*HBARC,
                       //		  RMS[0], RMS[1], R_WS);
                       point->a[count]=t;
                       point->b[count]=n_b;
                       point->c[count]=nmz;
                       point->d[count]=ia;
                       point->e[count]=N0_WS;
                       point->f[count]=Z0_WS;
                       point->g[count]=BEA_NUC*HBARC;
                       point->h[count]=BSA_NUC;
                       point->i[count]=BCA_NUC*HBARC;
                       point->j[count]=RMS[0];
                       point->k[count]=RMS[1];
                       point->l[count]=R_WS;
                           point->m[count]=ftime;


                       //		fprintf(myfile,"%e %e %i %i %e %e %e %e %e %e %e %e\n",point->a[count],point->b[count],point->c[count]
                       //		,point->d[count],point->e[count],point->f[count],point->g[count],point->h[count],point->i[count],point->j[count],point->k[count],
                       //		point->l[count])  ;

                       count++;




                       ic2 += 1;
                   }
                   }
                    else
                    {
                        ic = 0;
                    }*/
        }
    }


    /* last row
    nmz = ia = 0;
    A0_WS = N0_WS = Z0_WS = BEA_NUC = 0.;
    fprintf(FILE_OUT_EXTRA," %i %i %e %e %e %e\n",
      nmz, ia, A0_WS, N0_WS, Z0_WS, BEA_NUC*HBARC);
    fprintf(FILE_OUT_EXTRA,"# %i\n", ic2);
    */


    if(debug==1) fclose(FILE_PLOT2);
    if(debug==1)fclose(FILE_PLOT);
    //fclose(FILE_OUT_EXTRA);

    return point;
}
/*****************************************************************************/
int get_eos_table(int inse,int in_g,int in_e,int in_mu,int in_tau,int in_nu,
                  int *in_cl,int in_hyp,int in_tmes, int in_ton,
                  int *in_phase,int in_r, int in_tab,char *cnsfile)
{
    double t_min,t_max,t,n_b_min,n_b_max,y_q_min,y_q_max,n_b,y_q,d_y_q,
           f_n_b,m,dum,t_ref,y_q_ref,n_b_ref,f_t,inx_phase[N_PH];
    int n_n_b,i_n_b,n_y_q,i_y_q,i_y_q2,i_t,n_t,itmp,n_ibl,iwr,ip,iph,
        i_y_q_min,i_y_q_max,i_t_min,i_t_max,i_n_b_min,i_n_b_max,nnbp,ntp,
        in_cell,ibl,itmp1,itmp2,
        ir_y,ir_n,ir_t,idum,n_b_shift,
        n_para,i_para,ir;
    float dum_para[3];

    iwr = 0;
    FILE * TE;
    TE=fopen("XXXXXX.txt","a+");
    fclose(TE);
    IN_R = in_r;

    /* maximum number of blocks */
    n_ibl = 7;

    /* considered phases */
    in_phase[0] = 1;
    if (in_ton == 1)
    {
        for (iph=1; iph<N_PH; iph++)
        {
            in_phase[iph] = 0;
        }
    }

    /* only drop calculation stable? */
    in_cell  = 0;

    for (iph=1; iph<N_PH; iph++)
    {
        if (in_phase[iph] == 1) in_cell = 1;
    }

    /* reference mass = neutron mass [MeV] */
    m = M_NEU;
    m = 0.5*(M_NEU+M_PRO);

    /* read data file eos_table.in */
    if (in_tab == 1)
    {
        FILE_IN = fopen(changestring("eos_table",cnsfile,".in"),"r");
        fscanf(FILE_IN,"%i", &itmp);
        ibl = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_y_q_min = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_y_q_max = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_n_b_min = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_n_b_max = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_t_min = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_t_max = itmp;
    }
    else
    {
        FILE_IN = fopen(changestring("eos_tab",cnsfile,".in"),"r");
        fscanf(FILE_IN,"%i", &itmp);
        i_y_q_min = itmp;
        i_y_q_max = itmp;
        fscanf(FILE_IN,"%i", &itmp1);
        itmp = (itmp1-1)%50;
        ibl = (itmp1-itmp)/50+1;
        /*printf(" itmp ibl itmp1 %i %i %i\n", itmp1, ibl, itmp);*/
        i_n_b_min = itmp;
        i_n_b_max = itmp;
        fscanf(FILE_IN,"%i", &itmp);
        i_t_min = itmp;
        i_t_max = itmp;
    }
    fclose(FILE_IN);

    if (IN_R == 1)
    {
        for (i_n_b=0; i_n_b<DIM_N_B; i_n_b++)
        {
            for (i_y_q=0; i_y_q<DIM_Y_Q; i_y_q++)
            {
                for (iph=0; iph<N_PH; iph++)
                    R_IN[i_n_b][i_y_q][iph] = 0.;
            }
        }
        /* read data file eos_r1.in */
        FILE_IN = fopen(changestring("eos_r1",cnsfile,".in"),"r");
        do
        {
            fscanf(FILE_IN,"%i %i %f", &itmp1, &itmp2, &dum_para[0]);
            R_IN[itmp1][itmp2][1] = (double)dum_para[0];
            /*printf(" R1 %i %i %e\n", itmp1, itmp2, R_IN[itmp1][itmp2][1]);*/
        }
        while ((itmp1 > 0) && (itmp2 > 0));
        fclose(FILE_IN);
        /* read data file eos_r2.in */
        FILE_IN = fopen(changestring("eos_r2",cnsfile,".in"),"r");
        do
        {
            fscanf(FILE_IN,"%i %i %f", &itmp1, &itmp2, &dum_para[0]);
            R_IN[itmp1][itmp2][2] = (double)dum_para[0];
            /*printf(" R2 %i %i %e\n", itmp1, itmp2, R_IN[itmp1][itmp2][2]);*/
        }
        while ((itmp1 > 0) && (itmp2 > 0));
        fclose(FILE_IN);
        /* read data file eos_r3.in */
        FILE_IN = fopen(changestring("eos_r3",cnsfile,".in"),"r");
        do
        {
            fscanf(FILE_IN,"%i %i %f", &itmp1, &itmp2, &dum_para[0]);
            R_IN[itmp1][itmp2][3] = (double)dum_para[0];
            /*printf(" R3 %i %i %e\n", itmp1, itmp2, R_IN[itmp1][itmp2][3]);*/
        }
        while ((itmp1 > 0) && (itmp2 > 0));
        fclose(FILE_IN);
    }

    if ((iwr == 1) && (in_tab == 1))if(debug==1)fprintf(myfile,"\n block %i\n", ibl);

    /* charge fraction */
    if (i_y_q_min < 0) call_error(201);
    if (i_y_q_max > 100) call_error(202);
    /* Compose standard scaling */
    ir_y = 0;
    y_q_ref = 0.;
    d_y_q = 0.01;
    n_y_q = i_y_q_max-i_y_q_min+1;
    if ((iwr == 1) && (in_tab == 1))
    {
        y_q_min = d_y_q*(double)i_y_q_min;
        y_q_max = d_y_q*(double)i_y_q_max;
        if(debug==1)fprintf(myfile," range of charge fraction in calculation\n");
        if(debug==1)fprintf(myfile," Y_q_min = %f   Y_q_max = %f   # = %i\n",
                                y_q_min, y_q_max, n_y_q);
    }

    /* density */
    if (ibl > n_ibl) call_error(208);
    if (i_n_b_min < 0) call_error(205);
    if (i_n_b_max > 100) call_error(206);
    /* Compose standard */
    ir_n = 0;
    nnbp = 25;
    switch (ibl)
    {
    case 1:
    {
        n_b_ref = 1.e-12;
        n_b_shift = 1;
        break;
    }
    case 2:
    {
        n_b_ref = 1.e-10;
        n_b_shift = 2*nnbp+1;
        break;
    }
    case 3:
    {
        n_b_ref = 1.e-08;
        n_b_shift = 4*nnbp+1;
        break;
    }
    case 4:
    {
        n_b_ref = 1.e-06;
        n_b_shift = 6*nnbp+1;
        break;
    }
    case 5:
    {
        n_b_ref = 1.e-04;
        n_b_shift = 8*nnbp+1;
        break;
    }
    case 6:
    {
        n_b_ref = 1.e-02;
        n_b_shift = 10*nnbp+1;
        break;
    }
    case 7:
    {
        n_b_ref = 1.;
        n_b_shift = 12*nnbp+1;
        break;
    }
    default:
    {
        /* only for testing */
        n_b_ref = 0.4e-08; /* 0.4e-02 */

        n_b_ref = 0.01;

        nnbp = 25;
        n_b_shift = 0;

        n_b_ref = 1.e-10;
        nnbp = 2;
        n_b_shift = 0;
    }
    }
    f_n_b = pow(10.,(1./(double)nnbp));
    n_n_b = i_n_b_max-i_n_b_min+1;
    if ((iwr == 1) && (in_tab == 1))
    {
        n_b_min = n_b_ref*pow(f_n_b,(double)i_n_b_min);
        n_b_max = n_b_ref*pow(f_n_b,(double)i_n_b_max);
        if(debug==1)fprintf(myfile," range of baryon number density in calculation\n");
        if(debug==1)fprintf(myfile," n_b_min = %e fm^-3    n_b_max = %e fm^-3    # = %i\n",
                                n_b_min, n_b_max, n_n_b);
    }

    /* temperature */
    ir_t = 1;
    if (ir_t == 0)
    {
        /* Compose standard */
        t_ref = 0.1;
        ntp = 25;
        f_t = pow(10.,(1./(double)ntp));
        n_t = i_t_max-i_t_min+1;
        if (i_t_min < -1) call_error(203);
        if (i_t_max > 80) call_error(204);
        if (iwr == 1)
        {
            if (i_t_min == 0)
            {
                t_min = 0.;
            }
            else
            {
                t_min = t_ref*pow(f_t,(double)(i_t_min-1));
            }
            if (i_t_max == 0)
            {
                t_max = 0.;
            }
            else
            {
                t_max = t_ref*pow(f_t,(double)(i_t_max-1));
            }
        }
    }
    else
    {
        /* my grid */
        t_ref = 1.084202666;
        f_t = 0.09210340372;
        n_t = i_t_max-i_t_min+1;
        if (i_t_min < 0) call_error(203);
        if (i_t_max > 62) call_error(204);
        t_min = t_ref*sinh(f_t*(double)i_t_min);
        t_max = t_ref*sinh(f_t*(double)i_t_max);
    }
    if ((iwr == 1) && (in_tab == 1))
    {
        if(debug==1)fprintf(myfile," range of temperature in calculation\n");
        if(debug==1)fprintf(myfile," T_min = %e MeV   T_max = %e MeV   # = %i\n",
                                t_min, t_max, n_t);
    }

    /* standard tables for parameter-index correlation */
    get_parameter_tables(ir_y,ir_n,ir_t,d_y_q,f_n_b,t_ref,f_t,cnsfile);

    if(debug==1)FILE_OUT_THERMO = fopen(changestring("eos_thermo",cnsfile,".dat"),"w");
    if(debug==1)FILE_OUT_COMPO  = fopen(changestring("eos_compo",cnsfile,".dat"),"w");
    if(debug==1)FILE_OUT_COMPOH = fopen(changestring("eos_compo_hyp",cnsfile,".dat"),"w");
    FILE_OUT_EXTRA  = fopen(changestring("eos_extra",cnsfile,".dat"),"w");
    FILE_R1         = fopen(changestring("eos_r1",cnsfile,".dat"),"w");
    FILE_R2         = fopen(changestring("eos_r2",cnsfile,".dat"),"w");
    FILE_R3         = fopen(changestring("eos_r3",cnsfile,".dat"),"w");
    if(debug==1)FILE_PLOT       = fopen(changestring("plot",cnsfile,".dat"),"w");
    if(debug==1)FILE_PLOT2      = fopen(changestring("plot2",cnsfile,".dat"),"w");
    FILE_PARA       = fopen(changestring("para",cnsfile,".dat"),"r");

    /* non-existence of radius solutions before first run */
    for (ip=0; ip<N_PH; ip++) IEXR[ip] = 0;

    if (1 == 0)
    {
        /* read parameters from table */

        fscanf(FILE_PARA,"%i", &n_para);
        if(debug==1)fprintf(myfile," n_para %i\n", n_para);

        for (i_para=0; i_para<n_para; i_para++)
        {

            fscanf(FILE_PARA,"%i %e %e %e",
                   &idum, &dum_para[0], &dum_para[1], &dum_para[2]);
            t   = (double)dum_para[0];
            n_b = (double)dum_para[1]/CONV[0];
            y_q = (double)dum_para[2];
            if(debug==1)fprintf(myfile," %i %e %e %e\n", idum, t, n_b, y_q);

            /* considered phases */
            for (iph=0; iph<N_PH; iph++)
            {
                IN_PH[iph] = in_phase[iph];
                /*IHDB[iph] = in_phase[iph];*/
            }

            if ((y_q == 0.) || (y_q == 1.)) IN_PH[1] = IN_PH[2] = IN_PH[3] = 0;
            if (n_b > 0.15) IN_PH[1] = IN_PH[2] = IN_PH[3] = 0;
            if (n_b < 0.03) IN_PH[2] = IN_PH[3] = 0;

            for (ip=0; ip<N_PH; ip++)
            {
                IMFS[ip] = 0;
                /*IEXR[ip] = 0;*/
            }

            /***************************************************************/
            get_eos_local(t,n_b,y_q,inse,
                          in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                          in_cell,in_ton,0,0.,0.,cnsfile);
            /***************************************************************/

            if(debug==1)fprintf(FILE_PLOT," %i %e %e %e",
                                    idum, t, n_b, y_q);
            if(debug==1)fprintf(FILE_PLOT,
                                    " %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e\n",
                                    Y_PART[1], Y_PART[0],
                                    Y_PART[9], Y_PART[10], Y_PART[11], Y_PART[12]);


        }
        fclose(FILE_PARA);

    }
    else
    {

        /* hadronic charge fraction */
        for (i_y_q2=i_y_q_min; i_y_q2<(i_y_q_max+1); i_y_q2++)
        {
            /* hallo */

            if (i_y_q_min > 50)
            {
                i_y_q = i_y_q2;
                y_q = d_y_q*(double)i_y_q;

            }
            else
            {
                i_y_q = i_y_q_min-i_y_q2+i_y_q_max;
                y_q = d_y_q*(double)i_y_q;
            }

            /*d_y_q = 0.002;
            y_q = 0.03-d_y_q*(double)(i_y_q2-10);
            i_y_q = 10-i_y_q2;*/


            /*hallo*/
            /*y_q = 0.2055+0.0001*(double)i_y_q;*/

            /* considered phases */
            for (iph=0; iph<N_PH; iph++)
            {
                inx_phase[iph] = in_phase[iph];
                /*printf(" %i %i\n", iph, IN_PH[iph]);*/
                /*IHDB[iph] = in_phase[iph];*/
            }

            for (ip=0; ip<N_PH; ip++)
            {
                /*IMFS[ip] = 0;*/
                /* hallo */
                /*IEXR[ip] = 0;*/
            }

            /* baryonic number density */
            /* fit of flow
            i_n_b_min = 0;
            i_n_b_max = 532;*/
            for (i_n_b=i_n_b_min; i_n_b<(i_n_b_max+1); i_n_b++)
            {
                n_b = n_b_ref*pow(f_n_b,(double)i_n_b);
                /* hallo */
                /*n_b = 0.0011+0.00002*(double)i_n_b;*/
                /*n_b = 1.+0.02*(double)i_n_b;*/
                /*n_b = 0.1+0.02*(double)i_n_b;*/
                /*n_b = 0.0000001*(double)i_n_b;*/
                /*n_b = 0.0005*(double)i_n_b;*/
                /*n_b = 0.008+0.00002*(double)i_n_b;*/
                /*n_b = 4.e-03;*/
                /* for Rafal */
                n_b = 0.02*(double)i_n_b+1.;

                /* considered phases */
                for (iph=0; iph<N_PH; iph++)
                {
                    IN_PH[iph] = inx_phase[iph];
                    /*IHDB[iph] = in_phase[iph];*/
                }

                if (IN_R == 1)
                {

                    itmp = i_n_b+n_b_shift;
                    itmp1 = itmp+1;
                    if (itmp1 > (DIM_N_B-1)) itmp1 = DIM_N_B-1;
                    itmp2 = itmp-1;
                    if (itmp2 < 0) itmp2 = 0;

                    /* hallo
                    if ((inx_phase[1] == 1) &&
                        ((R_IN[itmp][i_y_q][1] > 0.) ||
                         (R_IN[itmp2][i_y_q][1] > 0.))) {
                      IN_PH[1] = 1;
                    }
                    else {
                      IN_PH[1] = 0;
                    }
                    if ((inx_phase[2] == 1) &&
                        ((R_IN[itmp][i_y_q][2] > 0.) ||
                         (R_IN[itmp1][i_y_q][2] > 0.) ||
                         (R_IN[itmp2][i_y_q][2] > 0.))) {
                      IN_PH[2] = 1;
                    }
                    else {
                      IN_PH[2] = 0;
                    }
                    if ((inx_phase[3] == 1) &&
                        ((R_IN[itmp][i_y_q][3] > 0.) ||
                         (R_IN[itmp1][i_y_q][3] > 0.) ||
                         (R_IN[itmp2][i_y_q][3] > 0.))) {
                      IN_PH[3] = 1;
                    }
                    else {
                      IN_PH[3] = 0;
                    }
                    */
                }

                if ((y_q == 0.) || (y_q == 1.)) IN_PH[1] = IN_PH[2] = IN_PH[3] = 0;
                if (n_b > 0.15) IN_PH[1] = IN_PH[2] = IN_PH[3] = 0;
                /* hallo
                if (n_b < 0.03) IN_PH[2] = IN_PH[3] = 0;
                    if (n_b > 0.07) IN_PH[1] = 0;
                    if (n_b > 0.125) IN_PH[2] = IN_PH[3] = 0;*/

                /* temperature */
                for (i_t=i_t_min; i_t<(i_t_max+1); i_t++)
                {
                    if (ir_t == 0)
                    {
                        if (i_t == 0)
                        {
                            t = 0.;
                        }
                        else
                        {
                            t = t_ref*pow(f_t,(double)(i_t-1));
                        }
                    }
                    else
                    {
                        t = t_ref*sinh(f_t*(double)i_t);
                    }

                    /*if (i_t == i_t_min) {
                      for (ip=0; ip<N_PH; ip++) IEXR[ip] = 0;
                      }*/

                    /*
                    	if (i_t == 0) {
                     	  if (i_n_b != i_n_b_min) {
                    	    if (ISOL == 2) IHDB[1] = 0;
                    	      if (ISOL == 3) IHDB[1] = IHDB[2] = 0;
                    	    if (ISOL == 0) IHDB[1] = IHDB[2] = IHDB[3] = 0;
                    	    if (ISOL == 0) in_phase[1] = in_phase[2] = in_phase[3] = 0;
                    	  }
                    	}
                    */

                    /* hallo */
                    if (ibl < 8)
                    {
                        /*t = 5.+0.05*(double)i_t;*/
                        t = 0.5*(double)i_t;
                        /*t = 5.*(double)i_t;*/
                        /* for Rafal */
                        t = 10.*(double)i_t;
                    }
                    /*t = 0.;*/

                    if (t < 270.)
                    {

                        IDX_TAB[0] = i_y_q;
                        IDX_TAB[1] = i_n_b+n_b_shift;
                        IDX_TAB[2] = i_t;


                        for (ip=0; ip<N_PH; ip++)
                        {
                            IMFS[ip] = 0;
                            /* hallo */
                            /*if (n_b > 0.1) IEXR[ip] = 0;*/
                        }

                        if(debug==1)fprintf(myfile,"\n *** y_q n_b t %e %e %e   %i %i %i\n",
                                                y_q, n_b, t, i_y_q, i_n_b, i_t);

                        if (1 == 0)
                        {
                            /* beta equilibrium */

                            /***************************************************************/
                            y_q = get_eos_local_beta(t,n_b,inse,
                                                     in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                                                     in_cell,in_ton,0,0.,0.,cnsfile);
                            /***************************************************************/

                            /*get_results_ws_cell();*/

                        }
                        else
                        {
                            /* hallo 2011/11/17 */
                            /*n_b = n_b*10.;*/
                            /***************************************************************/
                            get_eos_local(t,n_b,y_q,inse,
                                          in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                                          in_cell,in_ton,0,0.,0.,cnsfile);
                            /***************************************************************/

                        }
                        /*printf(" ICV_RMF %i\n", ICV_RMF);*/
                        if (ICV_RMF == 1)
                        {
                            /* compose standard */

                            /* thermodynamical properties */

                            if (1 == 0)
                            {
                                if(debug==1)fprintf(FILE_OUT_THERMO," %3i %3i %3i",
                                                        i_t, (i_n_b+n_b_shift), i_y_q);
                                if(debug==1)fprintf(FILE_OUT_THERMO," %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                                                        (HBARC*PRES/N_B),
                                                        S_DENS/N_B,
                                                        (HBARC*MU_BB/m-1.),
                                                        (HBARC*(MU_QQ)/m),
                                                        (HBARC*MU_LL/m));
                            }
                            else   /* for beta-equilibrium EoS */
                            {
                                if(debug==1)fprintf(FILE_OUT_THERMO,
                                                        " %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                                                        N_B, X_PART[15],  /*N_B*/
                                                        (HBARC*E_DENS/N_B-m), /*F*/
                                                        HBARC*PRES/N_B,  /*not /N_B*/
                                                        S_DENS/N_B, /*(HBARC*MU_BB-m), */
                                                        HBARC*(MU_QQ),
                                                        HBARC*MU_LL);

                                /*if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e\n",
                                    N_B, HBARC*PRES,
                                    (HBARC*F_DENS/N_B-m),
                                    (HBARC*E_DENS/N_B-m),
                                    (2.*DENS[9][0][0]/N_B));*/

                                /* for Rafal */
                                if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e\n",
                                                        t, N_B, HBARC*MU_BB, HBARC*PRES, HBARC*F_DENS/N_B);


                                /*if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e %e %e\n",
                                    N_B, HBARC*PRES, (HBARC*F_DENS/N_B-m),
                                    (DENS[26][0][0]+DENS[26][0][1])/N_B,
                                    DENS[27][0][0]/N_B, DENS[27][0][1]/N_B,
                                    (DENS[30][0][0]+DENS[30][0][1])/N_B);*/
                                /*if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e %e %e\n",
                                    N_B, HBARC*PRES, (HBARC*F_DENS/N_B-m),
                                    DENS[28][0][0]/N_B, DENS[28][0][1]/N_B,
                                    DENS[29][0][0]/N_B, DENS[29][0][1]/N_B);*/

                            }

                            /*
                            tmpn = exp((MU_BB-PARTICLE[1].m)/TT);
                            tmpp = exp((MU_BB+MU_QQ-PARTICLE[0].m)/TT);

                            if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e %e %e %e %e %e %e\n",
                                N_B, HBARC*TT,
                                X_PART[1],
                                S_DENS*TT/(2.5*PRES
                            	       -(MU_BB-PARTICLE[1].m)
                            	       *2.*tmpn*sqrt(CUBE((PARTICLE[1].m*TT)/TPI))
                            	       -(MU_BB+MU_QQ-PARTICLE[0].m)
                            	       *2.*tmpp*sqrt(CUBE((PARTICLE[0].m*TT)/TPI))),
                                (HBARC*E_DENS/N_B-m),
                                (HBARC*F_DENS/N_B-m),
                                HBARC*PRES/N_B,
                                S_DENS/N_B,
                                (HBARC*MU_BB-m),
                                HBARC*MU_QQ,
                                tmpn);*/
                            /*HBARC*MU_LL);*/

                            /*if(debug==1)fprintf(FILE_PLOT2, " %e %e %e %e %e %e %e\n",
                                N_B, HBARC*TT,
                                XP0/TT, XP1/TT, XP2/TT, MU_BB*HBARC,
                                (MU_BB-PARTICLE[1].m+SS[1][0]-VV[1][0])*HBARC);*/

                            if(debug==1)fprintf(FILE_PLOT2, " %e %e %e %e %e %e %e\n",
                                                    N_B, HBARC*TT,
                                                    X_PART[9], X_PART[13], X_PART[14], X_PART[15], X_PART[16]);


                            /*
                            if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e %e %e %e\n",
                                N_B, HBARC*TT, X_PART[0], X_PART[1],
                                X_NUC[0], X_NUC[1], X_NUC[2], X_NUC[3]);
                            */

                            fprintf(FILE_OUT_EXTRA,
                                    " %14.7e %14.7e %3i %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                                    N_B, y_q, (i_n_b+n_b_shift), R_WS, AA_WS, A0_WS, Z0_WS,
                                    HBARC*(MU_QQ+MU_LL));

                            /* composition */

                            if (1 == 0)
                            {
                                if(debug==1)fprintf(FILE_OUT_COMPO," %3i %3i %3i",
                                                        i_t, (i_n_b+n_b_shift), i_y_q);
                                if(debug==1)fprintf(FILE_OUT_COMPO,
                                                        " %2i %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e\n",
                                                        ISOL, A0_WS, Z0_WS,
                                                        Y_H,
                                                        Y_PART[1],
                                                        Y_PART[0],
                                                        Y_PART[2],
                                                        Y_PART[3],
                                                        Y_PART[9],
                                                        Y_PART[10],
                                                        Y_PART[12],
                                                        Y_PART[19]);
                            }
                            else
                            {
                                if(debug==1)fprintf(FILE_OUT_COMPO,
                                                        " %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                                                        N_B, Y_PART[2], Y_PART[3], AA_WS, A0_WS, Z0_WS, R_WS);
                                if(debug==1)fprintf(FILE_OUT_COMPOH,
                                                        " %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                                                        N_B, Y_S, X_PART[20],
                                                        X_PART[21], X_PART[22], X_PART[23],
                                                        X_PART[24], X_PART[25]);
                            }

                            /*
                                (Y_PART[1]),
                                (Y_PART[0]),
                                Y_PART[2], Y_PART[3],
                                (Y_PART[9]+Y_PART[13]),
                                (Y_PART[10]+Y_PART[17]),
                                (Y_PART[11]+Y_PART[18]),
                                (Y_PART[12]+Y_PART[19]));
                            */

                            if (t == 0.)
                            {
                                if (IN_PH[1] == 1)
                                {
                                    fprintf(FILE_R1," %3i %3i %11.5e\n",
                                            (i_n_b+n_b_shift), i_y_q, SOL[1].r);
                                }
                                if (IN_PH[2] == 1)
                                {
                                    fprintf(FILE_R2," %3i %3i %11.5e\n",
                                            (i_n_b+n_b_shift), i_y_q, SOL[2].r);
                                }

                                if (IN_PH[3] == 1)
                                {
                                    fprintf(FILE_R3," %3i %3i %11.5e\n",
                                            (i_n_b+n_b_shift), i_y_q, SOL[3].r);
                                }
                                /*
                                if(debug==1)fprintf(FILE_PLOT2,
                                    " %3i %3i %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e \n",
                                    (i_n_b+n_b_shift), i_y_q,
                                    HBARC*SOL[0].f/N_B, HBARC*SOL[1].f/N_B,
                                    HBARC*SOL[2].f/N_B, HBARC*SOL[3].f/N_B,
                                    SOL[1].r, SOL[2].r, SOL[3].r);
                                */
                            }

                            /* hallo */
                            /*
                            n_n = n_p = 0.5*N_B;
                            get_be_cl3(3,n_p,n_n);
                            if(debug==1)fprintf(myfile," %f %f\n", PARTICLE[9].be*HBARC, PARTICLE[12].be/DEJ);
                            if(debug==1)fprintf(myfile," %e %e\n", DEJ*HBARC, DEJT);
                            */
                            /*
                            if(debug==1)fprintf(FILE_PLOT," %e %e %e %e %e\n",
                                N_B, TT*HBARC,
                                (X_PART[10]+X_PART[17]),
                                (X_PART[11]+X_PART[18]),
                                (X_PART[12]+X_PART[19])
                                );
                            */

                            /* scalar densities and mean fields, only nucleons */
                            /*
                            fprintf(FILE_OUT_EXTRA," %3i %3i %3i",
                                i_t, (i_n_b+n_b_shift), i_y_q);
                            fprintf(FILE_OUT_EXTRA,
                                " %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e \n",
                                DENSS[1][0][2], DENSS[0][0][2],
                                VV[0][0]*HBARC, VV[1][0]*HBARC,
                                SS[0][0]*HBARC, SS[1][0]*HBARC);
                            */

                            if (1 == 0)
                            {
                                if (t == 0.)
                                {
                                    /* hallo */
                                    if (1 == 0)
                                    {
                                        for (ir=0; ir<NR; ir++)
                                        {
                                            if(debug==1)fprintf(FILE_PLOT," %f %e %e %e %e\n",
                                                                    RVEC[1][ir], DENS[0][ir][2], DENS[1][ir][2],
                                                                    DENS[2][ir][2], DENS[3][ir][2]);
                                        }
                                    }
                                    else
                                    {
                                        if(debug==1)fprintf(FILE_PLOT," %e %e %i %e %e %e %e %e %e %e\n",
                                                                N_B, Y_Q, ISOL,
                                                                HBARC*F_DENS/N_B,
                                                                HBARC*PRES,
                                                                /*-S_DENS/N_B,*/
                                                                HBARC*(F_DENS+PRES),
                                                                HBARC*(MU_BB+Y_Q*(MU_QQ+MU_LL))*N_B,
                                                                R_WS, A0_WS, Z0_WS);
                                        /*HBARC*MU_BB, HBARC*MU_QQ, HBARC*MU_LL);*/
                                        /*R_WS, A0_WS, Z0_WS, X_H);*/
                                    }
                                }
                            }


                        }

                        /*
                        if (ISOL == 2) in_drop = 0;
                        if (ISOL == 3) in_drop = in_bubble = 0;
                        if (ISOL == 0) in_drop = in_bubble = in_hole = 0;
                        */

                        /*if (ISOL == 2) IHDB[1] = 0;*/
                        /*if (ISOL == 3) IHDB[1] = IHDB[2] = 0;*/
                        /*if (ISOL == 0) IHDB[1] = IHDB[2] = IHDB[3] = 0;*/
                        if (0 == 1)
                        {
                            if (t == 0.)
                            {
                                if (ISOL == 0)
                                {
                                    /*
                                    if (fabs(SOL[0].f-SOL[1].f) > 0.001*N_B)
                                    inx_phase[1] = 0;
                                         if (fabs(SOL[0].f-SOL[2].f) > 0.001*N_B)
                                    inx_phase[2] = 0;
                                         */
                                    /*inx_phase[1] = inx_phase[2] = 0;*/
                                    /*if (fabs(SOL[0].f-SOL[3].f) > 0.001*N_B)
                                    inx_phase[3] = 0;*/
                                }
                                if (ISOL == 2)
                                {
                                    if (fabs(SOL[2].f-SOL[1].f) > 0.001*N_B)
                                        inx_phase[1] = 0;
                                }
                                if (ISOL == 3)
                                {
                                    /*
                                    if (fabs(SOL[3].f-SOL[1].f) > 0.001*N_B)
                                    inx_phase[1] = 0;
                                         */
                                    inx_phase[1] = 0;
                                    if (fabs(SOL[3].f-SOL[2].f) > 0.001*N_B)
                                        inx_phase[2] = 0;
                                }
                                if (SOL[2].r > 40.) inx_phase[2] = 0;
                            }
                        }
                        /* t > 0 ? */


                    }
                }
            }
        }
    }

    /* last row */
    if (1 == 0)
    {
        idum = -9;
        dum = 0.;
        if(debug==1)fprintf(FILE_OUT_THERMO," %3i %3i %3i", idum, idum, idum);
        if(debug==1)fprintf(FILE_OUT_THERMO," %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                                dum, dum, dum, dum, dum);
        if(debug==1)fprintf(FILE_OUT_COMPO," %3i %3i %3i", idum, idum, idum);
        idum = 0;

        if(debug==1)fprintf(FILE_OUT_COMPO,
                                " %2i %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e\n",
                                idum, dum, dum, dum,
                                dum, dum, dum, dum,
                                dum, dum, dum, dum);
    }

    if(debug==1) fclose(FILE_PLOT2);
    if(debug==1) fclose(FILE_PLOT);
    fclose(FILE_R3);
    fclose(FILE_R2);
    fclose(FILE_R1);
    if(debug==1)fclose(FILE_OUT_THERMO);
    if(debug==1)fclose(FILE_OUT_COMPO);
    if(debug==1)fclose(FILE_OUT_COMPOH);
    fclose(FILE_OUT_EXTRA);

    return 0;
}
/*****************************************************************************/
double get_eos_local_beta(double t,double n_b,int inse,
                          int in_g,int in_e,int in_mu,int in_tau,int in_nu,
                          int *in_cl,int in_hyp, int in_tmes,
                          int in_cell,int in_ton,int inuc,double a,double z,char* cnsfile)
{
    double xp,fp,xm,fm,x0,f0,dx;
    int icnt,iwr;

    iwr = 1;

    dx = 1.;

    if (n_b > 0.028)
    {
        xp = 0.2;
    }
    else
    {
        xp = 0.1*(log(n_b)+20.);
        if (xp > 0.)
        {
            xp = 0.47-0.1*QUAD(xp);
        }
        else
        {
            xp = 0.47;
        }
    }
    xp = 0.45;

    /*xp = 0.030;*/
    if(debug==1)fprintf(myfile," *** y_p_max fixed to %f in calculation of beta equilibrium\n", xp);

    icnt = 0;
    /***************************************************************/
    get_eos_local(t,n_b,xp,inse,
                  in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                  in_cell,in_ton,0,0.,0.,cnsfile);
    /***************************************************************/
    fp = HBARC*(MU_QQ+MU_LL);
    if (iwr == 1)if(debug==1)fprintf(myfile,"* %i %e %e %e %e %e\n",
                                             icnt, xp, dx, fp, HBARC*MU_QQ, HBARC*MU_LL);
    if (iwr == 2) if(debug==1)fprintf(FILE_PLOT2,
                                              "* %i %e %e %e %e %e\n",
                                              icnt, xp, dx, fp, HBARC*MU_QQ, HBARC*MU_LL);

    xm = xp;
    fm = fp;

    if (1 == 1)
    {
        xm = 0.027;
        xm = 0.3;
        if(debug==1)fprintf(myfile," *** y_p_min fixed to %f\n", xm);
        icnt += 1;
        /***************************************************************/
        get_eos_local(t,n_b,xm,inse,
                      in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                      in_cell,in_ton,0,0.,0.,cnsfile);
        /***************************************************************/
        fm = HBARC*(MU_QQ+MU_LL);
        if (iwr == 1)if(debug==1)fprintf(myfile,"* %i %e %e %e %e %e\n",
                                                 icnt, xm, dx, fm, HBARC*MU_QQ, HBARC*MU_LL);
        if (iwr == 2)if(debug==1)fprintf(FILE_PLOT2,
                                                 "* %i %e %e %e %e %e\n",
                                                 icnt, xm, dx, fm, HBARC*MU_QQ, HBARC*MU_LL);
    }

    while (fm > 0.)
    {
        xp = xm;
        fp = fp;
        xm *= 0.5;
        dx = xp-xm;
        icnt += 1;
        /***************************************************************/
        get_eos_local(t,n_b,xm,inse,
                      in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                      in_cell,in_ton,0,0.,0.,cnsfile);
        /***************************************************************/
        fm = HBARC*(MU_QQ+MU_LL);
        if (iwr == 1)if(debug==1)fprintf(myfile,"* %i %e %e %e %e %e\n",
                                                 icnt, xm, dx, fm, HBARC*MU_QQ, HBARC*MU_LL);
        if (iwr == 2) if(debug==1)fprintf(FILE_PLOT2,
                                                  "* %i %e %e %e %e %e\n",
                                                  icnt, xm, dx, fm, HBARC*MU_QQ, HBARC*MU_LL);
    }


    if(debug==1)fprintf(myfile,"\n %e %e\n", xm, fm);
    if(debug==1)fprintf(myfile," %e %e\n\n", xp, fp);

    f0 = 1.;
    x0 = 1.;
    /*while ((fabs(f0) > 1.e-08) && (icnt < 100)) {*/
    do
    {
        dx = x0;
        /*
        if (fabs(xp-xm) < 0.001) {
          x0 = 0.5*(xm+xp);
        }
        else {
          x0 = (fm*xp-fp*xm)/(fm-fp);
        }
        */
        x0 = 0.5*(xm+xp);

        dx -= x0;
        icnt += 1;
        /***************************************************************/
        get_eos_local(t,n_b,x0,inse,
                      in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                      in_cell,in_ton,0,0.,0.,cnsfile);
        /***************************************************************/
        f0 = HBARC*(MU_QQ+MU_LL);
        if (iwr == 1)if(debug==1)fprintf(myfile,"* %i %e %e %e %e %e\n",
                                                 icnt, x0, dx, f0, HBARC*MU_QQ, HBARC*MU_LL);
        if (iwr == 2) if(debug==1)fprintf(FILE_PLOT2,
                                                  "* %i %e %e %e %e %e\n",
                                                  icnt, x0, dx, f0, HBARC*MU_QQ, HBARC*MU_LL);
        if (f0 > 0.)
        {
            xp = x0;
            fp = f0;
        }
        else
        {
            xm = x0;
            fm = f0;
        }
    }
    while ((fabs(dx) > 1.e-09) && (icnt < 100));
    /*} while (((fabs(f0) > 1.e-06) || (fabs(dx) > 1.e-09)) && (icnt < 100));*/

    /*if (icnt > 99) {
     if(debug==1)fprintf(myfile," error in subroutine get_eos_local_beta\n");
      exit(0);
      }*/

    return x0;
}
/*****************************************************************************/
int get_eos_local(double t,double n_b,double y_q,int inse,int in_g,int in_e,int in_mu,int in_tau,int in_nu,int *in_cl,int in_hyp,int in_tmes,int in_cell,int in_ton,int inuc,double a,double z,char *cnsfile)
{
    int ip;
    double n_b_switch;

    /* 2011/11/18 */
    n_b_switch = 0.5;

    /* parameters */
    TT  = t/HBARC;
    N_B = n_b;
    Y_Q = y_q;

    Y_S = 0.;

    /* considered particles */
    IN_PART[0] = IN_PART[1] = 1;
    for (ip=2; ip<N_PART; ip++) IN_PART[ip] = 0;
    IN_PART[2] = in_e;
    IN_PART[3] = in_mu;
    IN_PART[4] = in_tau;
    IN_PART[5] = in_g;
    if (in_nu == 1)
    {
        for (ip=0; ip<3; ip++) IN_PART[6+ip] = 1;
    }
    /* clusters and nuclei, not for zero temperature or for large densities */
    if ((t > 0.) && (n_b < n_b_switch))   /* attention: depends on statistics and T */
    {
        for (ip=0; ip<N_CL; ip++) IN_PART[9+ip] = in_cl[ip];
        IN_TON = in_ton;
    }
    else
    {
        for (ip=0; ip<N_CL; ip++) IN_PART[9+ip] = 0;
        IN_TON = 0;
    }
    /* no mixed clusters or nuclei for neutron matter*/
    if (y_q < 0.5e-08)
    {
        IN_TON = 0;
        for (ip=0; ip<N_CL; ip++)
        {
            if (PARTICLE[9+ip].z > 0)
                IN_PART[9+ip] = 0;
        }
        /* no leptons for neutron matter */
        IN_PART[2] = IN_PART[3] = IN_PART[4] = 0;
    }
    if ((1.-y_q) < 0.5e-08)
    {
        IN_TON = 0;
        for (ip=0; ip<N_CL; ip++)
        {
            if (PARTICLE[9+ip].n > 0)
                IN_PART[9+ip] = 0;
        }
    }
    IN_CL = 0;
    for (ip=0; ip<N_CL; ip++)
    {
        if (IN_PART[9+ip] == 1) IN_CL = 1;
    }

    /* hyperons */
    if (in_hyp == 1)
    {
        IN_HYP = 1;
        for (ip=0; ip<N_HYP; ip++) IN_PART[9+N_CL+ip] = 1;
    }
    else
    {
        IN_HYP = 0;
    }

    /* thermal mesons */
    if ((in_tmes == 1) && (TT > 0.))
    {
        for (ip=0; ip<N_TMES; ip++) IN_PART[9+N_CL+N_HYP+ip] = 1;
        IN_TMES = 1;
    }
    else
    {
        for (ip=0; ip<N_TMES; ip++) IN_PART[9+N_CL+N_HYP+ip] = 0;
        IN_TMES = 0;
    }

    init_be_cl();

    if (inuc == 1)
    {
        /* calculation of finite nucleus */
        IN_TON  = 0;
        IN_CL   = 0;
        IN_HYP  = 0;
        for (ip=0; ip<N_HYP; ip++) IN_PART[9+N_CL+ip] = 0;
        IN_TMES = 0;
        for (ip=0; ip<N_TMES; ip++) IN_PART[9+N_CL+N_HYP+ip] = 0;
        if((int)get_nucleus(a,z,cnsfile) == 999) return 999;
        /*printf(" ICV_RMF = %i\n", ICV_RMF);*/
        /*printf(" %f %f %f %f %f %f %f %f %f %f\n",
           a, A0_WS, N0_WS, Z0_WS, R_WS, HBARC*BEA_NUC, HBARC*BFA_NUC, BSA_NUC,
           RMS[0], RMS[1]);*/

        if(debug==1)fprintf(myfile," %f %f %f %e %e %e %e %e\n",
                                N0_WS, Z0_WS, A0_WS,
                                BEA_NUC*HBARC, BSA_NUC,
                                RMS[0], RMS[1], R_WS);

    }
    else
    {
        /* calculation of equation of state */
        /* cell calculation always with electrons and without nuclei */
        if (in_cell == 1)
        {
            IN_PART[2] = 1;
            IN_TON = 0;
            IN_HYP = 0;
            for (ip=0; ip<N_HYP; ip++) IN_PART[9+N_CL+ip] = 0;
            IN_TMES = 0;
            for (ip=0; ip<N_TMES; ip++) IN_PART[9+N_CL+N_HYP+ip] = 0;
        }
        if (inse == 1)
        {
            /* NSE calculation always homogeneous */
            in_cell = 0;
            /*for (icl=0; icl<N_CL; icl++) IN_C[icl] = 0;
            IN_P[0] = IN_P[8] = IN_P[3] = 0;
            IN_P[11] = 1;*/
            /*get_eos_composition_nse(n_b,y_q,t);*/
        }
        else
        {
            if((int)get_eos_composition(in_cell,cnsfile)==999)
                return 999;
            /*get_properties();*/
            /*get_fractions();*/
        }
    }


    /*if (y_q <= 0.) X_B[0] = 0.;
      if (y_q >= 1.) X_B[1] = 0.;*/

    return 0;
}
/*****************************************************************************/
int get_ton_g(int ip)
{
    double fac,gx,dgx;

    fac = 0.;

    gx  = fac*QUAD(TT*HBARC);
    dgx = fac*2.*TT*QUAD(HBARC);

    NUC[ip].g = NUC[ip].g0+gx;
    NUC[ip].dlngdlnt = TT*dgx/NUC[ip].g;

    return 0;
}
/*****************************************************************************/
int get_ton_mass_shift(int ip,double n_ref,double n_eff)
{
    double tmp,f,df,fac;

    /* density dependence */
    fac = 0.*(double)NUC[ip].a/HBARC;

    tmp = n_eff/n_ref;
    f = tmp*(1.+tmp);
    df = (1.+2.*tmp)/n_ref;

    DM_TON = fac*f;

    DDMDN_TON = fac*df;

    /* temperature dependence */
    fac = 0./HBARC;

    f = fac*QUAD(TT*HBARC);
    df = fac*2.*TT*QUAD(HBARC);

    DM_TON += f;
    DDMDT_TON = df;

    return 0;
}
/*****************************************************************************/
int get_ton_mass(double n_p,double n_n)
{
    int ip;
    double n_eff,n_ref;

    n_ref = 0.15;

    for (ip=0; ip<N_NUC; ip++)
    {
        n_eff = 2.*(n_n*(double)NUC[ip].n+n_p*(double)NUC[ip].z)/(double)NUC[ip].a;

        get_ton_mass_shift(ip,n_ref,n_eff);
        NUC[ip].m = NUC[ip].m0+DM_TON;
        NUC[ip].dmdo = DDMDN_TON*LA_OMEGA;
        NUC[ip].dmdr = DDMDN_TON*NUC[ip].b*LA_RHO;
        NUC[ip].dmdt = DDMDT_TON;

        get_ton_g(ip);
    }

    return 0;
}
/*****************************************************************************/
int get_nucleus(double a,double z,char *cnsfile)
{
    int iwr,iph; /*,ir;*/
    double r,tmp1,tmp2,tmp1p,tmp2p,tmp1pp,tmp2pp,n,tmp2ppp;

    iwr = 0;

    n = a-z;
    AA_NUC = a;
    if (iwr == 1)if(debug==1)fprintf(myfile," AA_NUC %e\n", AA_NUC);

    init_sol();
    /* homogeneous matter */
    iph = 0;
    if (iwr == 1)if(debug==1)fprintf(myfile," homogeneous matter");
    IMFS[iph] = 0;
    iexs_sol(iph);
    r = 0.;
    XCL = 0;
    if (N_B > 0.)
    {
        if((int)solve_rmf(iph,r) == 999) return 999;
    }
    else
    {
        ISOL = 0;
        ICV_RMF = 1;
        NR = 1;
        MU_BB = PARTICLE[1].m;
        MU_QQ = PARTICLE[0].m-PARTICLE[1].m;
        MU_LL = PARTICLE[2].m;
        R_WS = r;
        /*printf(" m_p m_n %e %e\n", PARTICLE[0].m*HBARC, PARTICLE[1].m*HBARC);*/
        E_DENS = F_DENS = z*(PARTICLE[0].m+PARTICLE[2].m)+n*PARTICLE[1].m;
        C_DENS = 0.;
    }
    save_sol(iph);
    if (N_B > 0.)
    {
        tmp1   = E_DENS/N_B;
        tmp1p  = F_DENS/N_B;
        tmp1pp = S_DENS/N_B;
    }
    else
    {
        tmp1   = E_DENS/(double)a;
        tmp1p  = F_DENS/(double)a;
        tmp1pp = S_DENS/(double)a;
    }
    if (iwr == 1)if(debug==1)fprintf(myfile," iph F/A %i %f\n", iph, HBARC*tmp1);
    /* nucleus */
    iph = 1;
    if (iwr == 1)if(debug==1)fprintf(myfile," nucleus");
    IMFS[iph] = 0;
    iexs_sol(iph);
    if (N_B > 0.)
    {
        r = pow(3.*(double)a/(FPI*N_B),0.333333333333);
    }
    else
    {
        r = 20.;
    }
    if (iwr == 1)if(debug==1)fprintf(myfile," with cell radius %e fm\n", r);
    XCL = 0;
    if (N_B > 0.)
    {
        /*fit_radius_nuc(iph,r);*/
        if((int)solve_rmf(iph,r) == 999) return 999;
    }
    else
    {
        if((int)solve_rmf(iph,r) == 999) return 999;
    }

    save_sol(iph);


    if (N_B > 0.)
    {
        tmp2    = E_DENS/N_B;
        tmp2p   = F_DENS/N_B;
        tmp2pp  = S_DENS/N_B;
        tmp2ppp = C_DENS/N_B;
    }
    else
    {
        tmp2    = E_DENS*V_WS/AA_WS;
        tmp2p   = F_DENS*V_WS/AA_WS;
        tmp2pp  = S_DENS*V_WS/AA_WS;
        tmp2ppp = C_DENS*V_WS/AA_WS;
    }

    /* cm correction */
    if ((A0_WS > 0.) && (N_B == 0.))
    {
        /*tmp2 += -30.75/(HBARC*pow(A0_WS,1.333333333333));*/
        /*tmp2 += -17.2/(HBARC*pow(A0_WS,1.2));*/
    }




    if (iwr == 1)if(debug==1)fprintf(myfile," iph F/A %i %f\n", iph, HBARC*tmp2);

    if (iwr == 1)if(debug==1)fprintf(myfile," A0 N0 Z0 BE %f %f %f %f\n",
                                             A0_WS, N0_WS, Z0_WS, HBARC*(tmp1-tmp2));

    /*
    for (ir=0; ir<NR; ir++) {
      if(debug==1)fprintf(FILE_PLOT," %e %e %e %e\n", RVEC[1][ir],
        DENS[0][ir][2], DENS[1][ir][2], (DENS[0][ir][2]+DENS[1][ir][2]));
    }
    for (ir=0; ir<NR; ir++) {
      if(debug==1)fprintf(FILE_PLOT2," %e %e %e %e %e\n", RVEC[1][ir],
        SS[0][ir]*HBARC, SS[1][ir]*HBARC, VV[0][ir]*HBARC, VV[1][ir]*HBARC);
    }
    */

    BEA_NUC = tmp1-tmp2;
    BFA_NUC = tmp1p-tmp2p;
    BSA_NUC = tmp1pp-tmp2pp;
    BCA_NUC = -tmp2ppp;


    get_results_ws_cell(cnsfile);

    return 0;
}
/*****************************************************************************/
int fit_radius_nuc(int iph, double r0)
{
    int iwr,ir,ic,ic2;
    double r[3],f[3],dr,a,b,c,xp,xm,f0,da,da2,tmp;

    iwr = 1;

    dr = 1.;
    da = 1.e-03;
    da2 = QUAD(da);
    ic = 0;
    ic2 = 0;
    /*r0 += dr;*/

    for (ir=0; ir<3; ir++) r[ir] = f[ir] = 0.;

    if((int)solve_rmf(iph,r0) == 999) return 999;
    r[0] = r0;
    tmp = A0_WS-AA_NUC;
    f[0] = QUAD(tmp);
    if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r[0], f[0], tmp, A0_WS);

    if ((fabs(f[0]) > da2) && (TT > 0.))
    {
        r0 += dr;
        if((int)solve_rmf(iph,r0) == 999) return 999;
        r[1] = r0;
        tmp = A0_WS-AA_NUC;
        f[1] = QUAD(tmp);
        if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r[1], f[1], tmp, A0_WS);

        /*
        while (f[1] > f[0]) {
          r[2] = r[1];
          f[2] = f[1];
          r0 = 0.5*(r[0]+r[2]);
          dr*= 0.5;
          solve_rmf(iph,r0);
          r[1] = r0;
          f[1] = QUAD(A0_WS-AA_NUC);
          if (iwr == 1)if(debug==1)fprintf(myfile," %f %e\n", r[1], f[1]);
        }
        */
        while (f[1] > f[0])
        {
            r[2] = r[1];
            f[2] = f[1];
            r[1] = r[0];
            f[1] = f[0];
            r0 -= 2.*dr;
            if (r0 < 0.) r0 = 0.5*r[1];
            if((int)solve_rmf(iph,r0) == 999) return 999;
            r[0] = r0;
            tmp = A0_WS-AA_NUC;
            f[0] = QUAD(tmp);
            if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r[0], f[0], tmp, A0_WS);
        }

        if (r[2] == 0.)
        {
            r0 += dr;
            if((int)solve_rmf(iph,r0) == 999) return 999;
            r[2] = r0;
            tmp = A0_WS-AA_NUC;
            f[2] = QUAD(tmp);
            if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r[2], f[2], tmp, A0_WS);
            while (f[2] < f[1])
            {
                r[0] = r[1];
                f[0] = f[1];
                r[1] = r[2];
                f[1] = f[2];
                r0 += dr;
                if((int)solve_rmf(iph,r0) == 999) return 999;
                r[2] = r0;
                tmp = A0_WS-AA_NUC;
                f[2] = QUAD(tmp);
                if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r[2], f[2], tmp, A0_WS);
            }
        }


        if (iwr == 1)
        {
            if(debug==1)fprintf(myfile,"\n %f %e\n", r[0], f[0]);
            if(debug==1)fprintf(myfile," %f %e\n", r[1], f[1]);
            if(debug==1)fprintf(myfile," %f %e\n\n", r[2], f[2]);
        }


        while (ic == 0)
        {
            a = f[1];
            xp = r[2]-r[1];
            xm = r[0]-r[1];
            c = ((xm*f[2]-xp*f[0])/(xp-xm)+a)/(xm*xp);
            b = (f[2]-a-c*QUAD(xp))/xp;
            /*
            if(debug==1)fprintf(myfile,"\n %e %e\n", r[0], (a+xm*(b+xm*c)));
            if(debug==1)fprintf(myfile," %e %e\n", r[1], a);
            if(debug==1)fprintf(myfile," %e %e\n", r[2], (a+xp*(b+xp*c)));
            if(debug==1)fprintf(myfile,"\n a b c %e %e %e\n", a, b, c);
            */
            ic = 1;

            if (a > da2)
            {
                if (c/a > 0.)
                {
                    r0 = r[1]-0.5*b/c;
                    if (r0 < r[0]) r0 = 0.5*(r[0]+r[1]);
                    if (r0 > r[2]) r0 = 0.5*(r[1]+r[2]);
                    if((int)solve_rmf(iph,r0) == 999) return 999;
                    tmp = A0_WS-AA_NUC;
                    f0 = QUAD(tmp);
                    if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r0, f0, tmp, A0_WS);
                    ic2 = 1;
                    if (f0 > da2)
                    {
                        ic = 0;
                    }
                    if (r0 < r[1])
                    {
                        if (f0 < f[1])
                        {
                            r[2] = r[1];
                            f[2] = f[1];
                            r[1] = r0;
                            f[1] = f0;
                        }
                        else
                        {
                            r[0] = r0;
                            f[0] = f0;
                        }
                    }
                    else
                    {
                        if (f0 < f[1])
                        {
                            r[0] = r[1];
                            f[0] = f[1];
                            r[1] = r0;
                            f[1] = f0;
                        }
                        else
                        {
                            r[2] = r0;
                            f[2] = f0;
                        }
                    }
                }
            }

        }
        if (ic2 == 0)
        {
            r0 = r[1];
            if((int)solve_rmf(iph,r0) == 999) return 999;
            tmp = A0_WS-AA_NUC;
            f[1] = QUAD(tmp);
            if (iwr == 1)if(debug==1)fprintf(myfile," %f %e %f %f\n", r[1], f[1], tmp, A0_WS);
        }
    }


    return 0;
}
/*****************************************************************************/
int get_eos_composition(int in_cell,char* cnsfile)
{
    double r;
    int iwr,iph,cr,ip,ir;

    iwr = 1;

    init_sol();
    /* counter for phases */
    iph = 0;
    if (IN_PH[iph] == 1)
    {
        if (iwr == 1)if(debug==1)fprintf(myfile," homogeneous calculation\n");
        /*IPHASE = 0;*/
        /*IN_P[8] = IHNB = 0;*/
        IMFS[iph] = 0;
        iexs_sol(iph);
        r = 0.;
        XCL = 0;
        if((int)solve_rmf(iph,r)==999)
            return 999;
        save_sol(iph);
        if(debug==1)fprintf(myfile," F/A %f   %e\n", HBARC*F_DENS/N_B, DENS[12][0][2]);
    }
    if (in_cell == 1)
    {
        if (iwr == 1)if(debug==1)fprintf(myfile," inhomogeneous calculation\n");
        for (iph=1; iph<N_PH; iph++)
        {
            /*printf(" *** %i %i\n", iph, IN_PH[iph]);*/
            if (IN_PH[iph] == 1)
            {
                if (iwr)
                {
                    switch (iph)
                    {
                    case 2:
                    {
                        if(debug==1)fprintf(myfile," bubble \n");
                        break;
                    }
                    case 3:
                    {
                        if(debug==1)fprintf(myfile," hole \n");
                        break;
                    }
                    default:
                    {
                        if(debug==1)fprintf(myfile," drop \n");
                    }
                    }
                }
                /*IPHASE = 1;*/
                IMFS[iph] = 0;
                iexs_sol(iph);
                r = R_WS;
                XCL = 0;
                cr = fit_radius2(iph,r);
                if((int)cr == 999) return 999;
                save_sol(iph);
            }
        }
    }

    iph = select_sol(N_PH);
    if (iph < 0) return 0;

    iexs_sol(iph);

    /*
    if (iwr) {
     if(debug==1)fprintf(myfile," mu_b mu_q mu_l %e %e %e\n",
       MU_BB*HBARC, MU_QQ*HBARC, MU_LL*HBARC);
     if(debug==1)fprintf(myfile,"p  %e %e %e\n", DENS[0][0][0], DENS[0][0][1], DENS[0][0][2]);
     if(debug==1)fprintf(myfile,"n  %e %e %e\n", DENS[1][0][0], DENS[1][0][1], DENS[1][0][2]);
     if(debug==1)fprintf(myfile,"nn %e %e %e\n", DENS[15][0][0], DENS[15][0][1], DENS[15][0][2]);
    }
    */

    r = R_WS;
    XCL = 1;
    /*printf(" %i %i\n", iph, IMFS[iph]);*/
    /*IMFS[iph] = 0;*/
    if((int)solve_rmf(iph,r)==999)
        return 999;

    if (iwr == 1)if(debug==1)fprintf(myfile," iph ISOL R_WS Z0_WS A0_WS F/A  %i %i %f   %f %f   %f\n",
                                             iph, ISOL, R_WS, Z0_WS, A0_WS, F_DENS*HBARC/N_B);

    /*printf("* densities %e %e %e\n", DENS[0][0][2], DENS[1][0][2], DENS_TON[0]);*/

    if (ISOL == 0)
    {
        R_WS = A0_WS = Z0_WS = N0_WS = 0.;
        AA_WS = NN_WS = ZZ_WS = 0.;
    }

    get_fractions();

    /* Coulomb correction of chemical potentials */
    if ((ISOL > 0) && (IREXT == 1)) get_coul_corr();

    /* folding of cluster densities, for output */
    if (IPHASE > 0)
    {
        for (ip=0; ip<N_PART; ip++)
        {
            if (IN_PART[ip] == 1)
            {
                if (PARTICLE[ip].rms > 0.)
                {
                    /* vector density */
                    for (ir=0; ir<NR; ir++)
                    {
                        DENS_FOLD[ir] = DENS[ip][ir][2];
                    }
                    fold_dens(ip);
                    for (ir=0; ir<NR; ir++)
                    {
                        DENS[ip][ir][2] = DENS_FOLD[ir];
                    }
                    /* scalar density */
                    for (ir=0; ir<NR; ir++)
                    {
                        DENS_FOLD[ir] = DENSS[ip][ir][2];
                    }
                    fold_dens(ip);
                    for (ir=0; ir<NR; ir++)
                    {
                        DENSS[ip][ir][2] = DENS_FOLD[ir];
                    }
                }
            }
        }
    }

    if (IPHASE > 0)
    {
        if (ISOL > 0) get_results_ws_cell(cnsfile);
        /*printf(" get_results_ws_cell called\n");*/
    }

    return 0;
}
/*****************************************************************************/
int init_sol(void)
{
    int iph;

    for (iph=0; iph<N_PH; iph++)
    {
        SOL[iph].type  = 0;
        SOL[iph].conv  = 0;
        SOL[iph].nr    = 1;
        SOL[iph].f     = 0.;
        if (IEXR[iph] == 0)
        {
            SOL[iph].mu_b = 0.;
            SOL[iph].mu_q = 0.;
            SOL[iph].mu_l = 0.;
            SOL[iph].r    = 0.;
        }
    }

    return 0;
}
/*****************************************************************************/
int iexs_sol(int iph)
{
    double tmp;
    int ir,ip,iq,itmp,imf;

    IPHASE = iph;

    if (0 == 1)
    {
        if (iph == 0)
        {
            /* homogeneous */
            IPHASE = 0;
        }
        else
        {
            /* inhomogeneous */
            IPHASE = iph;
        }
    }

    /*printf(" iph IEXR %i %i\n", iph, IEXR[iph]);*/

    if (IMFS[iph] == 1)
    {
        ISOL    = SOL[iph].type;
        ICV_RMF = SOL[iph].conv;
        NR      = SOL[iph].nr;
        MU_BB   = SOL[iph].mu_b;
        MU_QQ   = SOL[iph].mu_q;
        MU_LL   = SOL[iph].mu_l;

        for (ip=0; ip<N_PART; ip++)
        {
            if (IN_PART[ip] == 1)
            {
                for (ir=0; ir<SOL[iph].nr; ir++)
                {
                    for (iq=0; iq<3; iq++)
                    {
                        DENS[ip][ir][iq] = XDENS[iph][ip][ir][iq];
                        DENSS[ip][ir][iq] = XDENSS[iph][ip][ir][iq];
                        DENSC[ip][ir][iq] = XDENSC[iph][ip][ir][iq];
                        DENSCS[ip][ir][iq] = XDENSCS[iph][ip][ir][iq];
                    }
                }
            }
        }
        for (imf=0; imf<N_MF; imf++)
        {
            for (ir=0; ir<SOL[iph].nr; ir++)
            {
                MF[imf][ir] = XMF[iph][imf][ir];
            }
        }

    }
    else
    {

        if ((N_B > 1.e-06) || (IPHASE > 0))
        {
            MU_BB = PARTICLE[1].m;
            MU_QQ = 0.;
        }
        else
        {
            if ((Y_Q < 1.) && (TT > 0.))
            {
                tmp = sqrt(TPI/(TT*PARTICLE[1].m));
                MU_BB = TT*log(0.5*N_B*(1.-Y_Q)*CUBE(tmp))+PARTICLE[1].m;
            }
            else
            {
                MU_BB = PARTICLE[1].m;
            }
            if ((Y_Q > 0.) && (TT > 0.))
            {
                tmp = sqrt(TPI/(TT*PARTICLE[0].m));
                MU_QQ = TT*log(0.5*N_B*Y_Q*CUBE(tmp))+PARTICLE[0].m-MU_BB;
            }
            else
            {
                MU_QQ = 0.;
            }
        }
        MU_LL = 0.;

    }

    if (IEXR[iph] == 1)
    {
        R_WS    = SOL[iph].r;
    }
    else
    {
        switch (iph)
        {
        case 0:
        {
            R_WS = 0.;
            break;
        }
        case 2:   /* bubble */
        {
            R_WS = 27.; /*27.;*/
            break;
        }
        case 3:   /* hole */
        {
            R_WS = 39.5; /*20.;*/
            break;
        }
        default:   /* drop */
        {
            /*
            R_WS = pow(14.4/N_B,0.333333333333);
            if (R_WS < 14.) {
            R_WS = 14.;
            }
                 */
            if (N_B > 1.e-03)
            {
                if (N_B > 0.03)
                {
                    R_WS = 16.4786+N_B*(-200.532+N_B*3121.1);
                }
                else
                {
                    R_WS = 23.68+log(N_B)*(6.506+log(N_B)*0.9947);
                }
                if (Y_Q > 0.3)
                {
                    R_WS *= 1.+1.25*(0.5-Y_Q);
                }
                else
                {
                    R_WS = pow((3.*(30.3*(1.-2.*Y_Q)+226.*Y_Q)
                                /(FPI*N_B*Y_Q)),0.333333333333);
                }
            }
            else
            {
                /*R_WS = exp(0.949242-0.327923*log(N_B));
                  R_WS *= 1.+1.31*(0.5-Y_Q);*/
                if (Y_Q > 0.62)
                {
                    R_WS = pow(4.54/(N_B*(1.-Y_Q)),0.333333333333);
                }
                else
                {
                    R_WS = pow(7.4/(N_B*Y_Q),0.333333333333);
                }

                /*R_WS = 90.4+0.*(0.1-Y_Q);*/
                /*R_WS = 194.;*/
            }
            /*R_WS = 20.;*/
        }
        }

        if (IN_R == 1)
        {
            tmp = R_IN[IDX_TAB[1]][IDX_TAB[0]][iph];
            if (tmp > 0.) R_WS = tmp;
            /*printf(" R_WS 0 %f\n", R_WS);*/
            if (R_WS < 0.1)
            {
                itmp = IDX_TAB[1]-1;
                if (itmp < 0) itmp = 0;
                tmp = R_IN[itmp][IDX_TAB[0]][iph];
                if (tmp > 0.)
                {
                    R_WS = tmp;
                    if (iph == 1) R_WS += 0.75;
                }
                /*printf(" R_WS - %f\n", R_WS);*/
                if (R_WS < 0.1)
                {
                    itmp = IDX_TAB[1]+1;
                    if (itmp > DIM_N_B) itmp = DIM_N_B;
                    tmp = R_IN[itmp][IDX_TAB[0]][iph];
                    /*printf(" R_WS + %f\n", R_WS);*/
                    if (tmp > 0.)
                    {
                        R_WS = tmp;
                        if (iph == 3) R_WS += 3.;
                    }
                }
            }
        }
    }

    return 0;
}
/*****************************************************************************/
int save_sol(int iph)
{
    int ip,ir,iq,imf;

    SOL[iph].type  = ISOL;
    SOL[iph].conv  = ICV_RMF;
    SOL[iph].nr    = NR;
    SOL[iph].mu_b  = MU_BB;
    SOL[iph].mu_q  = MU_QQ;
    SOL[iph].mu_l  = MU_LL;
    SOL[iph].r     = R_WS;
    SOL[iph].f     = F_DENS;

    if (ICV_RMF == 1)
    {
        IEXR[iph] = IMFS[iph] = 1;
        for (ip=0; ip<N_PART; ip++)
        {
            if (IN_PART[ip] == 1)
            {
                for (ir=0; ir<SOL[iph].nr; ir++)
                {
                    for (iq=0; iq<3; iq++)
                    {
                        XDENS[iph][ip][ir][iq] = DENS[ip][ir][iq];
                        XDENSS[iph][ip][ir][iq] = DENSS[ip][ir][iq];
                        XDENSC[iph][ip][ir][iq] = DENSC[ip][ir][iq];
                        XDENSCS[iph][ip][ir][iq] = DENSCS[ip][ir][iq];
                    }
                }
            }
        }
        for (imf=0; imf<N_MF; imf++)
        {
            for (ir=0; ir<SOL[iph].nr; ir++)
            {
                XMF[iph][imf][ir] = MF[imf][ir];
            }
        }
    }
    else
    {
        IEXR[iph] = IMFS[iph] = 0;
    }

    return 0;
}
/*****************************************************************************/
int select_sol(int ip_max)
{
    int ip,ip_min,ipp,iwr;
    double f_min;

    iwr = 1;

    for (ip=0; ip<ip_max; ip++)
    {
        if ((iwr == 1) && (IN_PH[ip] == 1))
            if(debug==1)fprintf(myfile," iph type conv F/A   %i   %i %i   %f\n",
                                    ip, SOL[ip].type, SOL[ip].conv, HBARC*SOL[ip].f/N_B);
        if ((ip > 0) && (fabs(SOL[ip].f-SOL[0].f) < 5.e-06*N_B))
        {
            SOL[ip].conv = 0;
        }
        if (SOL[ip].conv == 1)
        {
            for (ipp=0; ipp<ip; ipp++)
            {
                if ((SOL[ip].type == SOL[ipp].type) && (SOL[ipp].conv == 1) &&
                        (SOL[ip].f > SOL[ipp].f))
                {
                    SOL[ip].conv = 0;
                }
            }
        }
    }

    ip_min = -1;
    f_min = 1.e+04;
    for (ip=0; ip<ip_max; ip++)
    {
        if (SOL[ip].conv == 1)
        {
            IEXR[ip] = 1;
            if (SOL[ip].f < f_min)
            {
                ip_min = ip;
                f_min = SOL[ip].f;
            }
        }
        else
        {
            IEXR[ip] = 0;
        }
    }

    if (ip_min < 0)
    {
        ip_min = 0;
        ICV_RMF = 0;
    }

    return ip_min;
}
/*****************************************************************************/
int get_fractions(void)
{
    double tmpa,tmpv[DIM_R],tmpn,tmpz;
    int ip,ir;

    /*printf(" AA_WS A0_WS %f %f\n", AA_WS, A0_WS);*/
    for (ip=0; ip<N_PART; ip++)
    {
        if (ip < 5)
        {
            tmpa = 1.;
        }
        else
        {
            tmpa = (double)PARTICLE[ip].a;
        }
        X_PART[ip] = Y_PART[ip] = 0.;
        if (IN_PART[ip] == 1)
        {
            /* wozu st > 0 ? */
            if (PARTICLE[ip].st > 0)
            {
                if (ISOL == 0)   /* homogeneous */
                {
                    Y_PART[ip] = DENS[ip][0][2]/N_B;
                    X_PART[ip] = tmpa*Y_PART[ip];
                }
                else   /* inhomogeneous */
                {
                    for (ir=0; ir<NR; ir++)
                    {
                        /* needed because of folding */
                        DENS[ip][ir][2] = DENS[ip][ir][0]-DENS[ip][ir][1];
                        tmpv[ir] = DENS[ip][ir][2];
                    }
                    Y_PART[ip] = FPI*r_integrate(1,2,tmpv)/AA_WS;
                    X_PART[ip] = tmpa*Y_PART[ip];
                }
            }
        }
    }

    /*
    if(debug==1)fprintf(myfile," Z0 Z Z/A %f %e %e\n", Z0_WS, X_PART[0]*AA_WS, Z0_WS/AA_WS);
    if(debug==1)fprintf(myfile," N0 N N/A %f %e %e\n", N0_WS, X_PART[1]*AA_WS, N0_WS/AA_WS);
    if(debug==1)fprintf(myfile," A0_WS %f\n", A0_WS);
    */

    /* strangeness fraction */
    Y_S = 0.;
    if (ISOL == 0)
    {
        for (ip=20; ip<N_PART; ip++)
        {
            Y_S += (double)PARTICLE[ip].s*DENS[ip][0][2]/N_B;
        }
    }

    X_H = Y_H = 0.;

    /* drop/bubble correction */
    if ((ISOL == 1) || (ISOL == 2))
    {
        X_PART[0] -= Z0_WS/AA_WS;
        X_PART[1] -= N0_WS/AA_WS;
        Y_PART[0] = X_PART[0];
        Y_PART[1] = X_PART[1];
        Y_H = 1./AA_WS;
        X_H = A0_WS*Y_H;
    }

    /* table of nuclei */
    if (IN_TON == 1)
    {
        tmpn = tmpz = tmpa = 0.;
        for (ip=0; ip<N_NUC; ip++)
        {
            Y_NUC[ip] = DENS_TON[ip]/N_B;
            X_NUC[ip] = Y_NUC[ip]*(double)NUC[ip].a;
            tmpn += DENS_TON[ip]*(double)NUC[ip].n;
            tmpz += DENS_TON[ip]*(double)NUC[ip].z;
            tmpa += DENS_TON[ip];
        }
        X_H = (tmpn+tmpz)/N_B;
        if (tmpa > 0.)
        {
            N0_WS = tmpn/tmpa;
            Z0_WS = tmpz/tmpa;
            A0_WS = N0_WS+Z0_WS;
            Y_H = X_H/A0_WS;
        }
    }

    /*
    if(debug==1)fprintf(myfile," N0 Z0 A0 %f %f %f\n", N0_WS, Z0_WS, A0_WS);
    if(debug==1)fprintf(myfile," xp xn xa %e %e %e %e\n",
     X_PART[0], X_PART[1], X_NUC[0],
     (X_PART[0]+X_PART[1]+X_NUC[0]));
    if(debug==1)fprintf(myfile," yp yn ya %e %e %e\n", Y_PART[0], Y_PART[1], Y_NUC[0]);
    */

    /*
    if(debug==1)fprintf(myfile," X_p Z/A %e %e\n", X_PART[0], Z0_WS/AA_WS);
    if(debug==1)fprintf(myfile," X_n N/A %e %e\n", X_PART[1], N0_WS/AA_WS);
    if(debug==1)fprintf(myfile," A0_WS Z0_WS N0_WS %e %e %e\n", A0_WS, Z0_WS, N0_WS);
    if(debug==1)fprintf(myfile," Y_h X_h %e %e\n", Y_H, X_H);
    */

    /* limitation ??? */
    /*
    for (ip=0; ip<N_PART; ip++) {
      if (IN_PART[ip] == 1) {
        if (fabs(X_PART[ip]) < 1.e-10) X_PART[ip] = 0.;
      }
      }
    */

    /* effective fractions and densities, redistribution of two-body continuum states
    if (1 == 1) {
      if (IN_P[8] == 1) {
        for (ir=0; ir<NR; ir++) {
    DENS_B[0][ir][0] += (DENS_CL[5][ir]+2.*DENS_CL[7][ir]);
    DENS_B[1][ir][0] += (DENS_CL[5][ir]+2.*DENS_CL[6][ir]);
    DENS_CL[0][ir] += DENS_CL[4][ir];
        }
      }
      else {
        X_B[0] += 0.5*X_CL[5]+X_CL[7];
        X_B[1] += 0.5*X_CL[5]+X_CL[6];
        if (0 == 1) {
        X_CL[0] += X_CL[4];
    X_CL[4] = X_CL[5] = X_CL[6] = X_CL[7] = 0.;
        }
        X_CL[5] = X_CL[6] = X_CL[7] = 0.;
      }
    }
    */


    return 0;
}
/*****************************************************************************/
int get_coul_corr(void)
{
    int ip,ir;
    double fac,tmp0,ni,ne,tmp1[DIM_R];

    fac = QUAD(G_G)*RHO_MF[0][NRP]*(F_COUL2-F_COUL1);

    /* charged baryons */
    tmp0 = 0;
    for (ir=0; ir<NR; ir++) tmp1[ir] = 0.;
    /* protons */
    ip = 0;
    tmp0 += DENS[0][NRP][2];
    for (ir=0; ir<NR; ir++)
    {
        tmp1[ir] += DENS[ip][ir][2];
    }
    /* other particles */
    for (ip=9; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            tmp0 += (double)PARTICLE[ip].z*DENS[ip][NRP][2];
            for (ir=0; ir<NR; ir++)
            {
                tmp1[ir] += (double)PARTICLE[ip].z*DENS[ip][ir][2];
            }
        }
    }
    ni = FPI*r_integrate(0,2,tmp1);
    ne = VX*tmp0;
    MU_QQ += ne*fac/(ni+ne);

    /* charged leptons */
    tmp0 = 0.;
    for (ir=0; ir<NR; ir++) tmp1[ir] = 0.;
    for (ip=2; ip<5; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            tmp0 += DENS[ip][NRP][2];
            for (ir=0; ir<NR; ir++)
            {
                tmp1[ir] += DENS[ip][ir][2];
            }
        }
    }
    ni = FPI*r_integrate(0,2,tmp1);
    ne = VX*tmp0;
    MU_LL -= ne*fac/(ni+ne);

    return 0;
}
/*****************************************************************************/
int get_results_ws_cell(char *cnsfile)
{
    int ir;
    FILE *file_dat;

    if(debug==1)file_dat = fopen(changestring("dens",cnsfile,".dat"),"w");
    if(debug==1)fprintf(file_dat,"#           proton       neutron      total         ");
    if(debug==1)fprintf(file_dat,"isospin      electron     muon\n");
    if(debug==1)fprintf(file_dat,"# radius    density      density      density       ");
    if(debug==1)fprintf(file_dat,"density      density      density\n");
    if(debug==1)fprintf(file_dat,"#  [fm]     [fm^-3]      [fm^-3]      [fm^-3]       ");
    if(debug==1)fprintf(file_dat,"[fm^-3]      [fm^-3]      [fm^-3]\n");
    for (ir=0; ir<NR; ir++)
    {
        if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e %e %e %e %e\n", RVEC[1][ir],
                                DENS_N[0][ir], DENS_N[1][ir], DENS_N[2][ir], DENS_N[3][ir],
                                DENS[2][ir][2], DENS[3][ir][2],
                                DENS[9][ir][2], DENS[10][ir][2], DENS[12][ir][2], DENS[19][ir][2]);
    }
    /*
    if (IREXT == 1) {
      ir = NRP;
      if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e %e %e %e %e\n", R1,
        DENS_N[0][ir], DENS_N[1][ir], DENS_N[2][ir], DENS_N[3][ir],
        DENS[2][ir][2], DENS[3][ir][2],
        DENS[9][ir][2], DENS[10][ir][2], DENS[11][ir][2], DENS[12][ir][2]);
      if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e %e %e %e %e\n", R2,
        DENS_N[0][ir], DENS_N[1][ir], DENS_N[2][ir], DENS_N[3][ir],
        DENS[2][ir][2], DENS[3][ir][2],
        DENS[9][ir][2], DENS[10][ir][2], DENS[11][ir][2], DENS[12][ir][2]);
    }
    */
    if(debug==1) fclose(file_dat);

    if(debug==1)file_dat = fopen(changestring("denss",cnsfile,".dat"),"w");
    if(debug==1)fprintf(file_dat,"#           proton       neutron      total         ");
    if(debug==1)fprintf(file_dat,"isospin      electron     muon\n");
    if(debug==1)fprintf(file_dat,"# radius    density      density      density       ");
    if(debug==1)fprintf(file_dat,"density      density      density\n");
    if(debug==1)fprintf(file_dat,"#  [fm]     [fm^-3]      [fm^-3]      [fm^-3]       ");
    if(debug==1)fprintf(file_dat,"[fm^-3]      [fm^-3]      [fm^-3]\n");
    for (ir=0; ir<NR; ir++)
    {
        if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e\n", RVEC[1][ir],
                                DENSS_N[0][ir], DENSS_N[1][ir], DENSS_N[2][ir], DENSS_N[3][ir],
                                DENSS[2][ir][2], DENSS[3][ir][2]);
    }
    if(debug==1) fclose(file_dat);

    if(debug==1) file_dat = fopen(changestring("sources",cnsfile,".dat"),"w");
    if(debug==1)fprintf(file_dat,"# radius    gamma        omega        sigma        ");
    if(debug==1)fprintf(file_dat,"rho\n");
    if(debug==1)fprintf(file_dat,"#  [fm]     [fm^-3]      [fm^-3]      [fm^-3]      ");
    if(debug==1)fprintf(file_dat,"[fm^-3]\n");
    for (ir=0; ir<NR; ir++)
    {
        if(debug==1)fprintf(file_dat," %f %e %e %e %e\n", RVEC[1][ir],
                                S_MF[0][ir], S_MF[1][ir], S_MF[2][ir], S_MF[3][ir]);
    }
    if(debug==1) fclose(file_dat);

    if(debug==1) file_dat = fopen(changestring("mf",cnsfile,".dat"),"w");
    if(debug==1)fprintf(file_dat,"# radius    gamma        omega        sigma        ");
    if(debug==1)fprintf(file_dat,"rho\n");
    if(debug==1)fprintf(file_dat,"#  [fm]     [MeV]        [MeV]        [MeV]        ");
    if(debug==1)fprintf(file_dat,"[MeV]\n");
    for (ir=0; ir<NR; ir++)
    {
        if(debug==1)fprintf(file_dat," %f %e %e %e %e\n", RVEC[1][ir],
                                MF[0][ir]*HBARC, MF[1][ir]*HBARC, MF[2][ir]*HBARC,
                                MF[3][ir]*HBARC);
    }
    /*
    if (IREXT == 1) {
      ir = NRP;
      if(debug==1)fprintf(file_dat," %f %e %e %e %e\n", R1,
        MF[0][ir]*HBARC, MF[1][ir]*HBARC, MF[2][ir]*HBARC,
        MF[3][ir]*HBARC);
      if(debug==1)fprintf(file_dat," %f %e %e %e %e\n", R2,
        MF[0][ir]*HBARC, MF[1][ir]*HBARC, MF[2][ir]*HBARC,
        MF[3][ir]*HBARC);
    }
    */
    if(debug==1) fclose(file_dat);

    if(debug==1) file_dat = fopen(changestring("se",cnsfile,".dat"),"w");
    if(debug==1)fprintf(file_dat,"#           proton       proton       neutron      ");
    if(debug==1)fprintf(file_dat,"neutron       electron     muon\n");
    if(debug==1)fprintf(file_dat,"# radius    vector       scalar       vector       ");
    if(debug==1)fprintf(file_dat,"scalar        vector       vector\n");
    if(debug==1)fprintf(file_dat,"#  [fm]     [MeV]        [MeV]        [MeV]        ");
    if(debug==1)fprintf(file_dat,"[MeV]         [MeV]        [MeV]\n");
    for (ir=0; ir<NR; ir++)
    {
        if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e\n", RVEC[1][ir],
                                VV[0][ir]*HBARC, SS[0][ir]*HBARC,
                                VV[1][ir]*HBARC, SS[1][ir]*HBARC,
                                VV[2][ir]*HBARC, VV[3][ir]*HBARC);
    }
    /*
    if (IREXT == 1) {
      ir = NRP;
      if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e\n", R1,
        VV[0][ir]*HBARC, SS[0][ir]*HBARC,
        VV[1][ir]*HBARC, SS[1][ir]*HBARC,
        VV[2][ir]*HBARC, VV[3][ir]*HBARC);
      if(debug==1)fprintf(file_dat," %f %e %e %e %e %e %e\n", R2,
        VV[0][ir]*HBARC, SS[0][ir]*HBARC,
        VV[1][ir]*HBARC, SS[1][ir]*HBARC,
        VV[2][ir]*HBARC, VV[3][ir]*HBARC);
    }
    */
    if(debug==1) fclose(file_dat);

    return 0;
}
/*****************************************************************************/
int solve_rmf(int iph,double r)
{
    double mu_b,mu_q;
    int iwr,ic,ic3,ip;

    iwr = 0;

    /*ICON_CL = -1;*/

    IREXT = 0;
    if (IPHASE > 0)
    {
        R_WS = r;
        init_discr_ws_cell(iph);
        /*init_fourier_bessel();*/
        /*init_filon();*/
        if (R2 > R1) IREXT = 1;
    }
    else
    {
        NR  = 1;
        R_WS = 0.;
        V_WS = 1.;
    }
    NRP = NR-1;
    NX = NRP/2;

    if (N_B > 0.)
    {
        AA_WS = N_B*V_WS;
    }
    else
    {
        AA_WS = AA_NUC;
    }
    ZZ_WS = Y_Q*AA_WS;
    NN_WS = AA_WS-ZZ_WS;

    /*printf("***** iph IMFS %i %i\n", iph, IMFS[iph]);*/

    if (IMFS[iph] == 0)
    {
        /* initialisation of densities only if homogeneous calculation or
           not from previous calculation */
        init_dens();
        init_self_energies();
        init_mf();
        /* for pseudo densities, 2011/12/07: not here */
        /*if (LA_OMEGA != 0.) {
          for (ir=0; ir<NR; ir++) {
        MF[1][ir] = (DENS[0][ir][2]+DENS[1][ir][2])/LA_OMEGA;
          }
        }

        if (LA_RHO != 0.) {
          for (ir=0; ir<NR; ir++) {
        MF[3][ir] = (DENS[0][ir][2]-DENS[1][ir][2])/LA_RHO;
          }
          }*/
        get_self_energies(0);

    }
    else
    {
        /*init_self_energies();*/
        /*init_mf();*/
        get_self_energies(1);
    }

    /*get_test();*/

    /* initialisation of chemical potentials */
    mu_b = MU_BB;
    mu_q = MU_QQ;

    if ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1))
        MU_LL = get_chemical_potential_lep();
    /*printf(" MU_LL %e\n", MU_LL*HBARC);*/

    ICTR[0] = ICTR[1] = ICTR[2] = ICTR[3] = 0;

    /* first fit of chemical potentials:
       no relaxation of densities, no clusters, no hyperons */
    if (IMFS[iph] == 0) XHYP = 0;
    get_chemical_potential_bar(mu_b,mu_q);
    mu_b = F_MU_B;
    mu_q = F_MU_Q;

    /*printf(" mu_b mu_q %e %e\n", mu_b*HBARC, mu_q*HBARC);
      exit(0);*/

    /* initialize acceleration */
    /*init_acc();*/

    ic = 0;
    if (iwr == 1)if(debug==1)fprintf(myfile," *** %i %f %f %f\n",
                                             ic, mu_b*HBARC, mu_q*HBARC, MU_LL*HBARC);


    ic3 = 0;
    /*************************/
    if((int)cycle_solve(iph,mu_b,mu_q)==999)
        return 999;
    /*************************/
    /*mu_b = F_MU_B;
      mu_q = F_MU_Q;*/
    ic3 += CI_RMF;

    MU_BB = F_MU_B;
    MU_QQ = F_MU_Q;

    if (Y_Q == 0.)
    {
        if (TT > 0.)
        {
            MU_LL = 0.;
            MU_QQ = VV[0][0]-MU_BB;
        }
        else
        {
            MU_LL = PARTICLE[2].m;
            MU_QQ = VV[0][0]+fabs(PARTICLE[0].m-SS[0][0])-MU_BB;
        }
    }

    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            if ((ip < 2) || (ip > 8))
            {
                PARTICLE[ip].mu = MU_BB*(double)PARTICLE[ip].a
                                  +MU_QQ*(double)PARTICLE[ip].z;
                /*printf(" %i %i %i %e\n",
                  ip, PARTICLE[ip].a, PARTICLE[ip].z, PARTICLE[ip].mu*HBARC);*/
            }
            else
            {
                PARTICLE[ip].mu = MU_LL;
            }
        }
        else
        {
            PARTICLE[ip].mu = 0.;
        }
    }

    get_properties();

    get_anz_nuc();


    /* folding of cluster densities, for output */
    /* nicht hier
    if (IPHASE > 0) {
     if(debug==1)fprintf(myfile," fold\n");
      for (ip=0; ip<N_PART; ip++) {
        if (IN_PART[ip] == 1) {
    if (PARTICLE[ip].rms > 0.) {
      for (ir=0; ir<NR; ir++) {
        DENS_FOLD[ir] = DENS[ip][ir][2];
      }
      fold_dens(ip);
      for (ir=0; ir<NR; ir++) {
        DENS[ip][ir][2] = DENS_FOLD[ir];
      }
      for (ir=0; ir<NR; ir++) {
        DENS_FOLD[ir] = DENSS[ip][ir][2];
      }
      fold_dens(ip);
      for (ir=0; ir<NR; ir++) {
        DENSS[ip][ir][2] = DENS_FOLD[ir];
      }
    }
        }
      }
    }
    */

    /*
    if ((IN_P[7] == 1) && (IN_P[8] == 1)) {
    for (i1=0; i1<N_CL; i1++) {
      if (IN_C[i1] == 1) {
    for (ir=0; ir<NR; ir++) {
      DENS_FOLD[ir] = DENS_CL[i1][ir];
    }
    fold_dens(i1);
    for (ir=0; ir<NR; ir++)
      DENS_CL[i1][ir] = DENS_FOLD[ir];
      }
    }
    }
    */

    return 0;
}
/*****************************************************************************/
int cycle_solve(int iph,double mu_b,double mu_q)
{
    int ic,ic_max,iwr,iacc;
    double mu_b_old,mu_q_old,dmu_b,dmu_q,tmp,err_max;

    ICOUNT = 11111;

    iwr = 0;

    ic_max = 250; /* 250 */
    /* hallo TM1 */
    if (NL == 1) ic_max = 400;

    XCL = 0;

    if (IPHASE > -1)
    {
        /* homogeneous calculation */
        err_max = 0.5e-09;
        METHOD = 1;
        if (IN_HYP == 1) METHOD = 0;
        iacc = 0;

        init_acc0();

        ic = 0;
        tmp = 1.;
        do
        {
            ic += 1;

            mu_b_old = mu_b;
            mu_q_old = mu_q;
            get_self_energies(ic);

            if ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1))
                MU_LL = get_chemical_potential_lep();

            if(get_chemical_potential_bar(mu_b,mu_q) == 999) return 999;
            mu_b = F_MU_B;
            mu_q = F_MU_Q;

            dmu_b = fabs(mu_b-mu_b_old);
            dmu_q = fabs(mu_q-mu_q_old);
            tmp = dmu_b+dmu_q;

            if (iwr == 1)
            {
                if(debug==1)fprintf(myfile," *** %i %f %f %f %e %i %i %i   %e   %i %i\n",
                                        ic, F_MU_B*HBARC, F_MU_Q*HBARC,
                                        MU_LL*HBARC, tmp*HBARC,
                                        CI_YQ, CI_NB, CI_NBYQ, DENS[15][0][2], iacc, METHOD);
            }

            if (ic < 8)
            {
                tmp = 1.;
                if((int)get_acc0()==999)
                    return 999;
            }
            else
            {
                XHYP = 1;
                XCL = 1; /* not tested yet */
                if (iacc == 0)
                {
                    init_acc();
                    iacc = 1;
                }
                else
                {
                    if((int)get_acc()==999)
                        return 999;
                    if ((IPHASE > 0) && (tmp < 1.e-08)) METHOD = 0;
                }
            }

        }
        while ((tmp > err_max) && (ic < ic_max));

    }
    else
    {
        /* inhomogeneous calculation */
        err_max = 0.5e-09;
        METHOD = 1;
    }

    /* convergence check */
    ICV_RMF = 1;

    if (ic > (ic_max-1))
    {
        if ((IPHASE > 0) && (TT == 0.)) err_max *= 10.;
        if ((tmp > 2.*err_max) || (CI_NB > (ic_max-1)) || (CI_YQ > 199))
        {
            ICV_RMF = 0;
        }
    }

    /*printf(" ICV_RMF3 = %i", ICV_RMF);
    if (ICV_RMF == 0) {
     if(debug==1)fprintf(myfile,"     %i %i %i\n", ic, CI_NB, CI_YQ);
      }*/

    CI_RMF = ic;


    return 0;








    /*
    if (IEXR[iph] == 0) {
      XCL = 0;
    }
    else {
      XCL = 1;
    }
    */
    XCL = 0;

    get_acc0();

    ic = 0;
    tmp = 1.;
    do
    {
        ic += 1;
        ICOUNT = 11111; /*ic;*/

        /*if (ic == 4) init_acc();*/

        /*printf("%i get_self_energies\n", ic);*/
        mu_b_old = mu_b;
        mu_q_old = mu_q;
        get_self_energies(ic);

        if ((IPHASE > 0) &&
                ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1)))
            MU_LL = get_chemical_potential_lep();

        get_chemical_potential_bar(mu_b,mu_q);
        mu_b = F_MU_B;
        mu_q = F_MU_Q;

        dmu_b = fabs(mu_b-mu_b_old);
        dmu_q = fabs(mu_q-mu_q_old);
        tmp = dmu_b+dmu_q;


        if (iwr == 1)
        {
            if(debug==1)fprintf(myfile," *** %i %f %f %f %e %i %i %i   %e   %i %i\n",
                                    ic, F_MU_B*HBARC, F_MU_Q*HBARC,
                                    MU_LL*HBARC, tmp*HBARC,
                                    CI_YQ, CI_NB, CI_NBYQ, DENS[15][0][2], XCL, METHOD);
        }

        if ((IPHASE > 0) && (tmp < 1.e-09)) METHOD = 0;

        if ((IN_TON == 0) && (XCL == 1))
        {
            get_acc(); /* immer !?! */
        }
        else
        {
            get_acc0();
        }

        if ((IN_CL == 1) && (XCL == 0) && (ic > 4))
        {
            if (tmp < 1.e-05)
            {
                XCL = 1;
                init_acc();
                if (iwr == 1)if(debug==1)fprintf(myfile," clusters switched on\n");
            }
        }
        else   /* ??? */
        {
            if ((IN_TON == 0) && (XCL == 0) && (ic > 4))
            {
                if (tmp < 1.e-02)   /* 1.e-04 */
                {
                    XCL = 1;
                    init_acc();
                }
            }
        }

        if ((IN_CL == 1) && (XCL != 1)) tmp = 1.;

        if (ic < 8) tmp = 1.;


        /*
        if ((IPHASE > 0) && (tmp < 3.*err_max)) {
          METHOD = 0;
        }
        */

        /*if (ic%10 == 0) init_acc();*/


    }
    while ((tmp > err_max) && (ic < ic_max));


    /*
    get_self_energies(ic);
    if ((IPHASE > 0) &&
        ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1)))
      MU_LL = get_chemical_potential_lep();
    get_chemical_potential_bar(mu_b,mu_q);
    mu_b = F_MU_B;
    mu_q = F_MU_Q;
    */


    /* convergence check */
    ICV_RMF = 1;

    if (ic > (ic_max-1))
    {
        if ((IPHASE > 0) && (TT == 0.)) err_max *= 10.;
        if ((tmp > 2.*err_max) || (CI_NB > (ic_max-1)) || (CI_YQ > 199))
        {
            ICV_RMF = 0;
        }
    }

    /*printf(" ICV_RMF3 = %i", ICV_RMF);
    if (ICV_RMF == 0) {
     if(debug==1)fprintf(myfile,"     %i %i %i\n", ic, CI_NB, CI_YQ);
      }*/

    CI_RMF = ic;

    return 0;
}
/*****************************************************************************/
int cycle_solve_old(int iph,double mu_b,double mu_q)
{
    int ic,ic_max,iwr;
    double mu_b_old,mu_q_old,dmu_b,dmu_q,tmp,err_max;



    if ((TT > 0.) && (IPHASE > 0))
    {
        iwr = 1;
    }
    else
    {
        iwr = 1;
    }

    METHOD = 1;
    /*if (IN_HYP == 1) METHOD = 0;*/
    /*if (IN_CL == 1) METHOD = 0;*/

    ic_max = 250;
    /* hallo TM1 */
    if (NL == 1) ic_max = 400;

    /*ic_max = 0;*/

    err_max = 0.5e-09;
    /*if (IPHASE > 0) err_max = 0.5e-08;*/

    /*
    if (IEXR[iph] == 0) {
      XCL = 0;
    }
    else {
      XCL = 1;
    }
    */
    XCL = 0;

    get_acc0();

    ic = 0;
    tmp = 1.;
    do
    {
        ic += 1;
        ICOUNT = 11111; /*ic;*/

        /*if (ic == 4) init_acc();*/

        /*printf("%i get_self_energies\n", ic);*/
        mu_b_old = mu_b;
        mu_q_old = mu_q;
        get_self_energies(ic);

        if ((IPHASE > 0) &&
                ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1)))
            MU_LL = get_chemical_potential_lep();

        get_chemical_potential_bar(mu_b,mu_q);
        mu_b = F_MU_B;
        mu_q = F_MU_Q;

        dmu_b = fabs(mu_b-mu_b_old);
        dmu_q = fabs(mu_q-mu_q_old);
        tmp = dmu_b+dmu_q;


        if (iwr == 1)
        {
            if(debug==1)fprintf(myfile," *** %i %f %f %f %e %i %i %i   %e   %i %i\n",
                                    ic, F_MU_B*HBARC, F_MU_Q*HBARC,
                                    MU_LL*HBARC, tmp*HBARC,
                                    CI_YQ, CI_NB, CI_NBYQ, DENS[15][0][2], XCL, METHOD);
        }

        if ((IPHASE > 0) && (tmp < 1.e-09)) METHOD = 0;

        if ((IN_TON == 0) && (XCL == 1))
        {
            get_acc(); /* immer !?! */
        }
        else
        {
            get_acc0();
        }

        if ((IN_CL == 1) && (XCL == 0) && (ic > 4))
        {
            if (tmp < 1.e-05)
            {
                XCL = 1;
                init_acc();
                if (iwr == 1)if(debug==1)fprintf(myfile," clusters switched on\n");
            }
        }
        else   /* ??? */
        {
            if ((IN_TON == 0) && (XCL == 0) && (ic > 4))
            {
                if (tmp < 1.e-02)   /* 1.e-04 */
                {
                    XCL = 1;
                    init_acc();
                }
            }
        }

        if ((IN_CL == 1) && (XCL != 1)) tmp = 1.;

        if (ic < 8) tmp = 1.;


        /*
        if ((IPHASE > 0) && (tmp < 3.*err_max)) {
          METHOD = 0;
        }
        */

        /*if (ic%10 == 0) init_acc();*/


    }
    while ((tmp > err_max) && (ic < ic_max));


    /*
    get_self_energies(ic);
    if ((IPHASE > 0) &&
        ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1)))
      MU_LL = get_chemical_potential_lep();
    get_chemical_potential_bar(mu_b,mu_q);
    mu_b = F_MU_B;
    mu_q = F_MU_Q;
    */


    /* convergence check */
    ICV_RMF = 1;

    if (ic > (ic_max-1))
    {
        if ((IPHASE > 0) && (TT == 0.)) err_max *= 10.;
        if ((tmp > 2.*err_max) || (CI_NB > (ic_max-1)) || (CI_YQ > 199))
        {
            ICV_RMF = 0;
        }
    }

    /*printf(" ICV_RMF3 = %i", ICV_RMF);
    if (ICV_RMF == 0) {
     if(debug==1)fprintf(myfile,"     %i %i %i\n", ic, CI_NB, CI_YQ);
      }*/

    CI_RMF = ic;

    return 0;
}
/*****************************************************************************/
int get_properties(void)
{
    int icl,ir,ip,ih;
    double tmp,tmp1[DIM_R],tmp2[DIM_R],tmp3[DIM_R],tmp4[DIM_R],tmp5[DIM_R],
           io_dens,ie_dens,is_dens,ig_dens,ip_dens,isc_dens,iec_dens,
           xo_dens,xe_dens,xs_dens,xg_dens,xp_dens,xsc_dens,xec_dens,
           n_n,n_p,tmpv[DIM_R],dens_cl[N_CL][DIM_R],
           denss_cl[N_CL][DIM_R],s_g,p_g;

    /*zzzzzzz*/
    /*get_test2();
      exit(0);*/

    /*if(debug==1)fprintf(FILE_PLOT2," %e %e %e %e %e\n", N_B,
      (PARTICLE[1].m-SS[1][0])*HBARC,
      (PARTICLE[15].m-SS[15][0])*HBARC,
      DENS[15][0][2], DENSS[15][0][2]);*/

    /*printf(" mu_p mu_n %e %e\n", PARTICLE[0].mu*HBARC, PARTICLE[1].mu*HBARC);*/

    io_dens = ie_dens = is_dens = ig_dens = ip_dens = 0.;
    xo_dens = xe_dens = xs_dens = xg_dens = xp_dens = 0.;
    isc_dens = xsc_dens = 0.;
    iec_dens = xec_dens = 0.;

    C_DENS = 0.;

    /* free nucleons */
    for (ip=0; ip<2; ip++)
    {
        get_prop_part(ip);
        ie_dens += F_IE;
        is_dens += F_IS;
        ig_dens += F_IG;
        ip_dens += F_IP;
        if (ip == 1) XP0 = F_IP;
        if (IREXT == 1)
        {
            xe_dens += F_XE;
            xs_dens += F_XS;
            xg_dens += F_XG;
            xp_dens += F_XP;
        }
    }

    /*printf(" 0 ie is ig %e %e %e\n", ie_dens, is_dens, ig_dens);
     if(debug==1)fprintf(myfile," 0 xe xs xg %e %e %e\n", xe_dens, xs_dens, xg_dens);*/

    /* electrons, muons, tauons */
    for (ip=2; ip<5; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            get_prop_part(ip);
            ie_dens += F_IE;
            is_dens += F_IS;
            ig_dens += F_IG;
            ip_dens += F_IP;
            if (IREXT == 1)
            {
                xe_dens += F_XE;
                xs_dens += F_XS;
                xg_dens += F_XG;
                xp_dens += F_XP;
            }
        }
    }

    /* clusters, without external contribution, condensate contribution?
       only Maxwell Boltzmann corrected so far*/
    Y_C = 0.;
    if (IN_CL == 1)
    {
        /* cluster fraction */
        for (ir=0; ir<NR; ir++) tmpv[ir] = 0.;
        for (icl=0; icl<N_CL; icl++)
        {
            ip = 9+icl;
            if (IN_PART[ip] == 1)
            {
                for (ir=0; ir<NR; ir++)
                {
                    dens_cl[icl][ir] = DENS[ip][ir][0]-DENS[ip][ir][1]
                                       +DENSC[ip][ir][0]-DENSC[ip][ir][1];
                    denss_cl[icl][ir] = DENSS[ip][ir][0]-DENSS[ip][ir][1]
                                        +DENSC[ip][ir][0]-DENSC[ip][ir][1];
                    tmpv[ir] += dens_cl[icl][ir]*(double)PARTICLE[ip].a;
                }
            }
            else
            {
                for (ir=0; ir<NR; ir++) dens_cl[icl][ir] = 0.;
            }
        }
        if (IPHASE > 0)
        {
            tmp = FPI*r_integrate(1,2,tmpv)/V_WS;
        }
        else
        {
            tmp = tmpv[0];
        }
        Y_C = tmp/N_B;
        if (Y_C != 0.)
        {
            for (icl=0; icl<N_CL; icl++)
            {
                ip = 9+icl;
                if (IN_PART[ip] == 1)
                {
                    get_prop_part(ip);
                    /*printf(" ip F_IE F_IS F_IG F_IP %i %e %e %e %e\n",
                      ip, F_IE, F_IS, F_IG, F_IP);*/
                    ie_dens += F_IE;
                    is_dens += F_IS;
                    ig_dens += F_IG;
                    ip_dens += F_IP;
                    if (ip == 15)
                    {
                        /*XP1 = F_IS-XP3;
                          XP2 = XP3;*/
                        XP1 = F_IP;
                    }
                    if (IREXT == 1)
                    {
                        xe_dens += F_XE;
                        xs_dens += F_XS;
                        xg_dens += F_XG;
                        xp_dens += F_XP;
                    }
                    /* entropy correction */
                    if (IPHASE > 0)
                    {
                        for (ir=0; ir<NR; ir++)
                        {
                            n_n = 0.5*(RHO_PS[1][ir]-RHO_PS[3][ir]);
                            n_p = 0.5*(RHO_PS[1][ir]+RHO_PS[3][ir]);
                            get_be_cl3(icl,n_p,n_n);
                            tmpv[ir] = (denss_cl[icl][ir]-XXX*dens_cl[icl][ir])
                                       *(DBE_CL[2][icl]+DBE_CL[3][icl]);
                        }
                        isc_dens -= FPI*r_integrate(0,2,tmpv);
                        if (IREXT == 1)
                        {
                            xsc_dens -= VX*tmpv[NRP];
                        }
                    }
                    else
                    {
                        n_n = 0.5*(RHO_PS[1][0]-RHO_PS[3][0]);
                        n_p = 0.5*(RHO_PS[1][0]+RHO_PS[3][0]);
                        get_be_cl3(icl,n_p,n_n);
                        tmp = (denss_cl[icl][0]-XXX*dens_cl[icl][0]) /*2011/11/10*/
                              *(DBE_CL[2][icl]+DBE_CL[3][icl]);
                        isc_dens -= tmp;
                    }
                    /* energy shift correction ??? */
                    if (IPHASE > 0)
                    {
                        for (ir=0; ir<NR; ir++)
                        {
                            n_n = 0.5*(RHO_PS[1][ir]-RHO_PS[3][ir]);
                            n_p = 0.5*(RHO_PS[1][ir]+RHO_PS[3][ir]);
                            get_be_cl3(icl,n_p,n_n);
                            tmpv[ir] = dens_cl[icl][ir]*BE_CL[icl];
                        }
                        iec_dens += FPI*r_integrate(0,2,tmpv);
                        if (IREXT == 1)
                        {
                            xec_dens += VX*tmpv[NRP];
                        }
                    }
                    else
                    {
                        n_n = 0.5*(RHO_PS[1][0]-RHO_PS[3][0]);
                        n_p = 0.5*(RHO_PS[1][0]+RHO_PS[3][0]);
                        get_be_cl3(icl,n_p,n_n);
                        tmp = dens_cl[icl][0]*BE_CL[icl];
                        iec_dens += tmp;
                    }
                }
            }
        }
    }

    /* hyperons */
    if (IN_HYP == 1)
    {
        for (ih=0; ih<N_HYP; ih++)
        {
            ip = 9+N_CL+ih;
            if (IN_PART[ip] == 1)
            {
                get_prop_part(ip);
                ie_dens += F_IE;
                is_dens += F_IS;
                ig_dens += F_IG;
                ip_dens += F_IP;
                if (IREXT == 1)
                {
                    xe_dens += F_XE;
                    xs_dens += F_XS;
                    xg_dens += F_XG;
                    xp_dens += F_XP;
                }
            }
        }
    }

    /* thermal mesons */
    if (IN_TMES == 1)
    {
        for (ih=0; ih<N_TMES; ih++)
        {
            ip = 9+N_CL+N_HYP+ih;
            if (IN_PART[ip] == 1)
            {
                get_prop_part(ip);
                ie_dens += F_IE;
                is_dens += F_IS;
                ig_dens += F_IG;
                ip_dens += F_IP;
                if (IREXT == 1)
                {
                    xe_dens += F_XE;
                    xs_dens += F_XS;
                    xg_dens += F_XG;
                    xp_dens += F_XP;
                }
            }
        }
    }

    /* interacting meson fields */
    if (IPHASE > 0)
    {
        for (ir=0; ir<NR; ir++)
        {
            /*tmp = DENS_N[2][ir];*/
            tmp = RHO_TOT[ir];
            get_cpl(tmp,0);
            /* free meson contribution, with partial integration and field eqs.,
            sources include thermal and condensate contributions */
            tmp1[ir] = MF[2][ir]*S_MF[2][ir]
                       +MF[4][ir]*S_MF[4][ir]
                       +MF[5][ir]*S_MF[5][ir]
                       -MF[1][ir]*S_MF[1][ir]
                       -MF[3][ir]*S_MF[3][ir];
            /* interaction contribution, only thermal part, no condensate! */
            /* xxxxxx RHO_MF */
            tmp2[ir] = CPL[1][0]*MF[1][ir]*RHO_MF[1][ir]
                       +CPL[3][0]*MF[3][ir]*RHO_MF[3][ir]
                       +CPL[5][0]*MF[5][ir]*RHO_MF[5][ir];
        }
        tmp1[0] = FPI*r_integrate(0,2,tmp1);
        tmp2[0] = FPI*r_integrate(0,2,tmp2);
        /* surface contributions */
        /*
          nrp = NRP;
          nrq = NR-2;
          tmp = FPI*RVEC[2][nrp]/(RVEC[1][nrp]-RVEC[1][nrq]);
          tmp1[0] -= MF[1][nrp]*(MF[1][nrp]-MF[1][nrq])*tmp;
          tmp2[0] += MF[2][nrp]*(MF[2][nrp]-MF[2][nrq])*tmp;
          tmp1[0] -= MF[3][nrp]*(MF[3][nrp]-MF[3][nrq])*tmp;
        */
        /*tmp1[0] /= V_WS;
          tmp2[0] /= V_WS;*/
        ie_dens += 0.5*tmp1[0]+tmp2[0];
        if (IREXT == 1)
        {
            /* free meson contribution */
            xe_dens += VX*(0.5*tmp1[NRP]+tmp2[NRP]);
        }
    }
    else
    {
        /*tmp = DENS_N[2][0];*/
        tmp = RHO_TOT[0];
        get_cpl(tmp,0);
        tmp1[0] = MESON[1].m2*QUAD(MF[1][0]);
        tmp2[0] = MESON[2].m2*QUAD(MF[2][0]);
        tmp3[0] = MESON[3].m2*QUAD(MF[3][0]);
        tmp4[0] = MESON[4].m2*QUAD(MF[4][0]);
        tmp5[0] = MESON[5].m2*QUAD(MF[5][0]);
        tmp = 0.5*(tmp2[0]+tmp4[0]-tmp1[0]-tmp3[0]-tmp5[0]);
        ie_dens += tmp;
        ip_dens -= tmp;
        XP2 = -tmp;

        if (NL == 1)
        {
            tmp1[0] = NL_GO*MF[1][0]*(RHO_MF[1][0]+RHOC_MF[1][0]);
            tmp3[0] = NL_GR*MF[3][0]*(RHO_MF[3][0]+RHOC_MF[3][0]);
            tmp5[0] = 0.;
        }
        else
        {
            tmp1[0] = CPL[1][0]*MF[1][0]*(RHO_MF[1][0]+RHOC_MF[1][0]);
            tmp3[0] = CPL[3][0]*MF[3][0]*(RHO_MF[3][0]+RHOC_MF[3][0]);
            tmp5[0] = CPL[5][0]*MF[5][0]*(RHO_MF[5][0]+RHOC_MF[5][0]);
        }
        ie_dens += tmp1[0]+tmp3[0]+tmp5[0];

        /* rearrangement contribution from nucleons */

        tmp1[0] = CPL[1][1]*MF[1][0]*(RHO_MF[1][0]+RHOC_MF[1][0])
                  -CPL[2][1]*MF[2][0]*(RHO_MF[2][0]+RHOC_MF[2][0])
                  +CPL[3][1]*MF[3][0]*(RHO_MF[3][0]+RHOC_MF[3][0])
                  -CPL[4][1]*MF[4][0]*(RHO_MF[4][0]+RHOC_MF[4][0]);
        ip_dens += tmp1[0]*DENS_N[2][0];
        XP2 += tmp1[0]*DENS_N[2][0];
        XP3 = 0.;

        if (NL == 1)
        {
            tmp2[0] = -NL_G2*CUBE(MF[2][0])/3.+NL_G3*QUAR(MF[2][0])/4.
                      -NL_C3*QUAR(MF[1][0])/4.;
            ie_dens += tmp2[0];
        }
    }

    /* Coulomb field */
    if (IPHASE > 0)
    {
        /* interior volume contribution */
        for (ir=0; ir<NR; ir++) tmp1[ir] = MF[0][ir]*RHO_MF[0][ir];
        tmp1[0] = G_G*FPI*r_integrate(0,2,tmp1);
        /* surface contribution */
        if (IREXT == 1)
        {
            tmp = QUAD(G_G*RHO_MF[0][NRP])*VX*F_COUL1;
            tmp1[0] -= tmp;
        }
        ie_dens += 0.5*tmp1[0];
        C_DENS = 0.5*tmp1[0];
        /* exterior contribution to Coulomb energy */
        if (IREXT == 1)
        {
            tmp = 0.5*VX*QUAD(G_G*RHO_MF[0][NRP])*(F_COUL2+F_COUL1);
            xe_dens += tmp;
            C_DENS = C_DENS+tmp;
        }
    }

    /* Coulomb correction */
    if (IREXT == 1)
    {
        DMU= VX*QUAD(G_G*RHO_MF[0][NRP])*(F_COUL2-F_COUL1);
        xg_dens += DMU;
    }
    else
    {
        DMU = 0.;
    }

    /* entropy correction from clusters */
    if (IN_CL == 1)
    {
        ie_dens += TT*isc_dens;
        is_dens += isc_dens;
        if (IREXT == 1)
        {
            xe_dens += TT*xsc_dens;
            xs_dens += xsc_dens;
        }
    }

    /* energy shift correction from clusters */
    /*if (IN_CL == 1) {
      ie_dens += iec_dens;
      if (IREXT == 1) {
        xe_dens += xec_dens;
      }
      }*/

    /* grand thermodynamical potential */
    io_dens = ie_dens-TT*is_dens-ig_dens;
    /*io_dens = - ip_dens;*/
    if ((IPHASE > 0) && (IREXT == 1))
    {
        xo_dens = xe_dens-TT*xs_dens-xg_dens;
    }

    /* thermodynamical properties */
    O_DENS = io_dens;
    S_DENS = is_dens;
    G_DENS = ig_dens;
    if (IPHASE > 0)
    {
        if (IREXT == 1)
        {
            O_DENS += xo_dens;
            S_DENS += xs_dens;
            G_DENS += xg_dens;
        }
        O_DENS /= V_WS;
        S_DENS /= V_WS;
        G_DENS /= V_WS;
        C_DENS /= V_WS;
    }

    /* table of nuclei */
    if (IN_TON == 1)
    {
        get_prop_ton();
        O_DENS += O_TON;
        S_DENS += S_TON;
        G_DENS += G_TON;
    }

    /* photons */
    if (IN_PART[5] == 1)
    {
        s_g = 4.*PI2*CUBE(TT)/45.;
        p_g = 0.25*TT*s_g;
        O_DENS -= p_g;
        S_DENS += s_g;
    }

    PRES = -O_DENS;
    /*PRES = ip_dens;*/
    /*XP3 = S_DENS;*/
    F_DENS = G_DENS+O_DENS;
    E_DENS = F_DENS+TT*S_DENS;

    /* energy modification 2011/11/18 */
    /*printf(" %e\n", (MM[15][0]-2.*PARTICLE[1].m)*HBARC*dens_cl[6][0]/N_B);
      E_DENS -= (MM[15][0]-2.*PARTICLE[1].m)*dens_cl[6][0];*/

    return 0;
}
/*****************************************************************************/
int get_prop_ton(void)
{
    int ip,iq;

    O_TON = S_TON = G_TON = 0.;
    for (ip=0; ip<NP_NUC; ip++)
    {
        iq = ISORT[ip];
        /*printf(" iq ip mu m VV SS dens %i %i %e %e %e %e %e\n",
           iq, ip, NUC[iq].mu*HBARC, NUC[iq].m*HBARC,
           VV_NUC[iq]*HBARC, SS_NUC[iq]*HBARC,
           exp((NUC[iq].mu-NUC[iq].ms-VV_NUC[iq])/TT));
          if(debug==1)fprintf(myfile," %e %e\n", DENS_TON[iq], DENSS_TON[iq]);*/
        O_TON += DENS_TON[iq];
        G_TON += NUC[iq].mu*DENS_TON[iq];
        S_TON += ((2.5+NUC[iq].dlngdlnt
                   -(NUC[iq].mu-VV_NUC[iq]-NUC[iq].ms)/TT)*DENS_TON[iq]
                  -NUC[iq].dmdt*DENSS_TON[iq]);
    }
    O_TON *= -TT;

    /*printf(" O G S %e %e %e\n", O_TON*HBARC, G_TON*HBARC, S_TON);
     if(debug==1)fprintf(myfile," E/A %e\n", (O_TON+G_TON+S_TON*TT)*HBARC/N_B);*/

    return 0;
}
/*****************************************************************************/
int get_anz_nuc(void)
{
    double max[3],min[3],tmpv[DIM_R],tmpv2[DIM_R],tmp[2],tmp2[2];
    int in,ir,ir_max[3],ir_min[3];

    RMS[0] = RMS[1] = 0.;

    if (IPHASE > 0)
    {

        for (ir=0; ir<NR; ir++)
        {
            DENS_N[0][ir] = DENS[0][ir][0]-DENS[0][ir][1];
            DENS_N[1][ir] = DENS[1][ir][0]-DENS[1][ir][1];
            DENS_N[2][ir] = DENS_N[0][ir]+DENS_N[1][ir];
        }

        for (in=0; in<3; in++)
        {
            /* maximum */
            ir_max[in] = 0;
            max[in] = 0.;
            for (ir=0; ir<NR; ir++)
            {
                if (DENS_N[in][ir] > max[in])
                {
                    ir_max[in] = ir;
                    max[in] = DENS_N[in][ir];
                }
            }
            /* minimum */
            ir_min[in] = ir_max[in];
            min[in] = max[in];
            for (ir=0; ir<NR; ir++)
            {
                if (DENS_N[in][ir] < min[in])
                {
                    ir_min[in] = ir;
                    min[in] = DENS_N[in][ir];
                }
            }
            /*
            if(debug==1)fprintf(myfile," in ir_min r min %i %3i %e %e\n",
            in, ir_min[in], RVEC[1][ir_min[in]], min[in]);
            if(debug==1)fprintf(myfile," in ir_max r max %i %i %e %e\n",
            in, ir_max[in], RVEC[1][ir_max[in]], max[in]);
                 */
        }
        if (ir_min[2] > ir_max[2])
        {
            /* drop */
            ISOL = 1;
            for (in=0; in<2; in++)
            {
                for (ir=0; ir<ir_min[in]; ir++)
                {
                    tmpv[ir] = DENS_N[in][ir]-min[in];
                    tmpv2[ir] = tmpv[ir]*RVEC[2][ir];
                }
                for (ir=ir_min[in]; ir<NR; ir++)
                {
                    tmpv[ir] = 0.;
                    tmpv2[ir] = 0.;
                }
                /*
                for (ir=0; ir<NR; ir++) {
                  if(debug==1)fprintf(FILE_PLOT," %f %e %e\n", RVEC[1][ir], tmpv[ir], tmpv2[ir]);
                }
                exit(0);
                */
                tmp[in] = FPI*r_integrate(0,2,tmpv);
                tmp2[in] = FPI*r_integrate(0,2,tmpv2);
                if (tmp[in] > 0.) RMS[in] = sqrt(tmp2[in]/tmp[in]);
                /* no external contribution */
            }
        }
        else
        {
            /*printf(" max %f %f\n", max[2], DENS_N[2][NRP]);*/
            if (DENS_N[2][NRP] < 0.75*max[2])
            {
                /* bubble */
                ISOL = 2;
                for (in=0; in<2; in++)
                {
                    if (ir_min[in] < ir_max[in])
                    {
                        for (ir=0; ir<ir_min[in]; ir++)
                        {
                            tmpv[ir] = 0.;
                        }
                        for (ir=ir_min[in]; ir<NR; ir++)
                        {
                            tmpv[ir] = DENS_N[in][ir]-min[in];
                        }
                    }
                    else
                    {
                        for (ir=0; ir<ir_min[in]; ir++)
                        {
                            tmpv[ir] = DENS_N[in][ir]-min[in];
                        }
                        for (ir=ir_min[in]; ir<NR; ir++)
                        {
                            tmpv[ir] = 0.;
                        }
                    }
                    tmp[in] = FPI*r_integrate(0,2,tmpv);
                    /* no external contribution */
                }
            }
            else
            {
                /* hole */
                ISOL = 3;
                for (in=0; in<2; in++)
                {
                    for (ir=0; ir<ir_max[in]; ir++)
                    {
                        tmpv[ir] = DENS_N[in][ir]-max[in];
                    }
                    for (ir=ir_max[in]; ir<NR; ir++)
                    {
                        tmpv[ir] = 0.;
                    }
                    tmp[in] = FPI*r_integrate(0,2,tmpv);
                    /* no external contribution */
                }
            }
        }
        Z0_WS = tmp[0];
        N0_WS = tmp[1];
    }
    else
    {
        /* homogeneous */
        ISOL = 0;
        Z0_WS = N0_WS = 0.;
    }
    A0_WS = Z0_WS+N0_WS;

    return 0;
}
/*****************************************************************************/
int get_test2()
{
    int ix;
    double mu,tmp,tmps;

    if(debug==1)fprintf(myfile,"hallo test2\n");
    IPHASE = 0;
    VV[2][0] = 0.;
    SS[2][0] = 0.;
    TT = 75./HBARC;
    PARTICLE[2].st = 1;
    PARTICLE[2].m = 500./HBARC;
    for (ix=-100; ix<101; ix++)
    {
        mu = PARTICLE[2].mu = 6.*(double)ix/HBARC;
        get_prop_part(2);
        tmp = F_IE;
        tmps = F_IS;
        if(debug==1)fprintf(FILE_PLOT," %e %e %e\n",
                                mu*HBARC, tmp, tmps);
    }

    return 0;
}
/*****************************************************************************/
int get_prop_part(int ip)
{
    double tmp1[DIM_R],tmp2[DIM_R],tmp3[DIM_R],tmp4[DIM_R],mu,d;
    int ir;

    F_ST = PARTICLE[ip].st;
    F_G  = PARTICLE[ip].g;
    F_M  = PARTICLE[ip].m;
    mu   = PARTICLE[ip].mu;

    /*if (ip>1)if(debug==1)fprintf(myfile," ip T g %i %e %e\n", ip, TT*HBARC, F_G);*/

    if (IPHASE > 0)
    {
        for (ir=0; ir<NR; ir++)
        {
            F_S = SS[ip][ir];
            F_V = VV[ip][ir];
            get_prop(mu);
            tmp1[ir] = F_ENE;
            tmp2[ir] = F_ENT;
            /* unfolded densities, particle - antiparticle */
            tmp3[ir] = DENS[ip][ir][0]-DENS[ip][ir][1];
            tmp4[ir] = F_PRE;
        }
        F_IE = FPI*r_integrate(0,2,tmp1);
        F_IS = FPI*r_integrate(0,2,tmp2);
        F_IG = mu*FPI*r_integrate(0,2,tmp3);
        F_IP = FPI*r_integrate(0,2,tmp4);
        if (IREXT == 1)
        {
            F_XE = VX*tmp1[NRP];
            F_XS = VX*tmp2[NRP];
            F_XG = mu*VX*tmp3[NRP];
            F_XP = VX*tmp4[NRP];
        }
    }
    else
    {
        F_S = SS[ip][0];
        F_V = VV[ip][0];
        get_prop(mu);
        F_IE = F_ENE;
        F_IS = F_ENT;
        F_IG = mu*(DENS[ip][0][0]-DENS[ip][0][1]);
        F_IP = F_PRE;
        F_XE = F_XS = F_XG = F_XP = 0.;
        /* entropy and energy correction
           because of temperature dependent degeneracy factor */
        if ((TT > 0.) && (PARTICLE[ip].tdgdt != 0.))
        {
            d = PARTICLE[ip].tdgdt*F_IP;
            F_IE += d;
            F_IS += d/TT;
            XP3 = d/TT;
        }
        /*
        if (ip > 1) {
        if(debug==1)fprintf(myfile," ip e e %i %e %e\n", ip, F_IE, (TT*F_IS-F_IP+F_IG-F_V*DENS[ip][0][2]));
        if(debug==1)fprintf(myfile," F_G s %e %e\n", F_G, (F_IE+F_IP-F_IG+F_V*DENS[ip][0][2])/TT);
        if(debug==1)fprintf(myfile," s ds %e %e\n", F_IS, PARTICLE[ip].tdgdt*F_IP/TT);
        }
        */
    }

    return 0;
}
/*****************************************************************************/
int get_prop(double mu)
{
    int ipa;
    double fac,ml,kf,tmp1,tmp2,ene,ent,pre,tmpn,tmpns,mup;

    F_MD = F_M-F_S;
    F_MD2 = QUAD(F_MD);
    fac = F_G/TPI2;

    F_ENE = F_ENT = F_PRE = 0.;

    if (TT > 0.)
    {
        switch(F_ST)
        {
        case -1:   /* Bose-Einstein, particle = antiparticle */
        {
            mup = mu;
            if (mu > (F_V+fabs(F_MD))) mup = F_V+fabs(F_MD);
            IBX = ipa = 0;
            get_simint_bounds_b(mup,ipa);
            F_ENE += fac*get_simint(10,1.e-08,14,F_KM,F_KP,mup,ipa);
            F_ENT -= fac*get_simint(10,1.e-08,16,F_KM,F_KP,mup,ipa);
            F_PRE += fac*get_simint(10,1.e-08,15,F_KM,F_KP,mup,ipa)/3.;
            break;
        }
        case 0:   /* Bose-Einstein, particle != antiparticle */
        {
            /* limits on chemical potential */
            mup = mu;
            if (mu > (F_V+fabs(F_MD))) mup = F_V+fabs(F_MD);
            if (mu < (F_V-fabs(F_MD))) mup = F_V-fabs(F_MD);
            for (ipa=0; ipa<2; ipa++)
            {
                IBX = ipa;
                get_simint_bounds_b(mup,ipa);
                F_ENE += fac*get_simint(10,1.e-08,14,F_KM,F_KP,mup,ipa);
                F_ENT -= fac*get_simint(10,1.e-08,16,F_KM,F_KP,mup,ipa);
                F_PRE += fac*get_simint(10,1.e-08,15,F_KM,F_KP,mup,ipa)/3.;
            }
            break;
        }
        case 2:   /* Maxwell-Boltzmann, non-relativistic */
        {
            /*
            tmp1 = sqrt(fabs(F_MD)*TT/TPI);
            tmp1 = CUBE(tmp1);
            tmp2 = tmp1*F_G*exp((mu-fabs(F_MD)-F_V)/TT);
            tmpn = tmp2;
            tmpns = tmp2*(1.-1.5*TT/fabs(F_MD));
            F_ENE = 3.*TT*tmpn+fabs(F_MD)*tmpns;
            F_ENT = 4.*tmpn+(fabs(F_MD)*tmpns-mu*tmpn)/TT;
            tmp2 = tmp1*F_G*exp((F_V-fabs(F_MD)-mu)/TT);
            tmpn = tmp2;
            tmpns = tmp2*(1.-1.5*TT/fabs(F_MD));
            F_ENE += 3.*TT*tmpn+fabs(F_MD)*tmpns;
            F_ENT += 4.*tmpn+(fabs(F_MD)*tmpns+mu*tmpn)/TT;
            */
            tmp1 = sqrt(fabs(F_MD)*TT/TPI);
            /*tmp1 = sqrt(fabs(F_M)*TT/TPI);*/
            tmp1 = CUBE(tmp1);
            tmpn = tmp1*F_G*exp((mu-fabs(F_MD)-F_V)/TT);
            F_ENE = (fabs(F_MD)+1.5*TT)*tmpn;
            F_ENT = tmpn*(2.5+(fabs(F_MD)-mu+F_V)/TT);
            F_PRE = tmpn*TT;
            /*printf(" (mu-m+F_S-F_V) mu F_G m S V F_PRE %e %e %e %e %e %e\n",
               (mu-fabs(F_MD)-F_V)*HBARC, F_G, F_M*HBARC,
               F_S*HBARC, F_V*HBARC, F_PRE*HBARC);*/
            /* hallo antiparticles
            tmpn = tmp1*F_G*exp((F_V-fabs(F_MD)-mu)/TT);
            F_ENE += (fabs(F_MD)+1.5*TT)*tmpn;
            F_ENT += tmpn*(2.5+(fabs(F_MD)+mu-F_V)/TT);
            F_PRE += tmpn*TT;
            */
            break;
        }
        case 3:   /* Maxwell-Boltzmann, relativistic */
        {
            tmp1 = sqrt(fabs(F_MD)*TT/TPI);
            tmp1 = F_G*CUBE(tmp1);
            tmp2 = tmp1*exp((mu-fabs(F_MD)-F_V)/TT);
            tmpn = tmp2*get_k2_red(fabs(F_MD)/TT);
            tmpns = tmp2*get_k1_red(fabs(F_MD)/TT);
            F_ENE = 3.*TT*tmpn+fabs(F_MD)*tmpns;
            F_ENT = 4.*tmpn+(fabs(F_MD)*tmpns-(mu-F_V)*tmpn)/TT;
            F_PRE = tmpn*TT;
            tmp2 = tmp1*exp((F_V-fabs(F_MD)-mu)/TT);
            tmpn = tmp2*get_k2_red(fabs(F_MD)/TT);
            tmpns = tmp2*get_k1_red(fabs(F_MD)/TT);
            F_ENE += 3.*TT*tmpn+fabs(F_MD)*tmpns;
            F_ENT += 4.*tmpn+(fabs(F_MD)*tmpns+(mu-F_V)*tmpn)/TT;
            F_PRE += tmpn*TT;
            break;
        }
        default:   /* Fermi-Dirac */
        {
            for (ipa=0; ipa<2; ipa++)
            {
                get_simint_bounds_f(mu,ipa);
                if (F_KP > F_KM)
                {
                    ene = get_simint(10,1.e-08,4,F_KM,F_KP,mu,ipa);
                    ent = get_simint(10,1.e-08,6,F_KM,F_KP,mu,ipa);
                    pre = get_simint(10,1.e-08,5,F_KM,F_KP,mu,ipa)/3.;
                }
                else
                {
                    ene = ent = pre = 0.;
                }
                if (F_KM > 0.)
                {
                    tmp1 = CUBE(F_KM)/3.;
                    ml = sqrt(F_MD2+QUAD(F_KM));
                    tmp1 *= ml;
                    ene += 0.75*tmp1;
                    pre += 0.25*tmp1;
                    tmp1 = 0.125*F_MD2*(F_KM*ml-F_MD2*log((F_KM+ml)/fabs(F_MD)));
                    ene += tmp1;
                    pre -= tmp1;
                }
                F_ENE += ene;
                F_ENT += ent;
                F_PRE += pre;
            }
            F_ENE *=  fac;
            F_ENT *= -fac;
            F_PRE *=  fac;
        }
        }
    }
    else
    {
        if (F_ST == 1)   /* Fermi-Dirac */
        {
            ml = fabs(mu-F_V);
            tmp1 = QUAD(ml)-F_MD2;
            if (tmp1 > 0.)
            {
                kf = sqrt(tmp1);
                tmp2 = ml*kf*tmp1/3.;
                F_ENE = 0.75*tmp2;
                F_PRE = 0.25*tmp2;
                tmp2 = 0.125*F_MD2*(kf*ml-F_MD2*log((kf+ml)/fabs(F_MD)));
                F_ENE += tmp2;
                F_PRE -= tmp2;
                F_ENE *= fac;
                F_PRE *= fac;
            }
        }
    }

    return 0;
}
/*****************************************************************************/
int get_acc(void)
{
    int ip,iq,ir,idx[REC_ACC],ik,in,itmp;
    double relax,tmpf,tmpm,dmax;

    for (ip=0; ip<REC_ACC; ip++)
    {
        INDEX[ip] = idx[ip] = (ip+M_ACC)%REC_ACC;
    }

    /*
    for (ip=0; ip<REC_ACC; ip++) {
      if ((idx[ip] < 0) || (idx[ip] > (REC_ACC+1))) {
       if(debug==1)fprintf(myfile," ip idx %i %i\n", idx[ip]);
        call_error(3333);
      }
    }
    */

    M_ACC += 1;

    iq = 0;
    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {

            if ((IPHASE > 0) && (PARTICLE[ip].rms > 0.))
            {
                /* vector density */
                for (ir=0; ir<NR; ir++)
                {
                    DENS_FOLD[ir] = DENS[ip][ir][2];
                }
                fold_dens(ip);
                for (ir=0; ir<NR; ir++)
                {
                    DENS[ip][ir][2] = DENS_FOLD[ir];
                }
                /* scalar density */
                for (ir=0; ir<NR; ir++)
                {
                    DENS_FOLD[ir] = DENSS[ip][ir][2];
                }
                fold_dens(ip);
                for (ir=0; ir<NR; ir++)
                {
                    DENSS[ip][ir][2] = DENS_FOLD[ir];
                }
            }

            for (ir=0; ir<NR; ir++)
            {
                VMB[iq][idx[0]] = DENS[ip][ir][2];
                iq += 1;
            }
            for (ir=0; ir<NR; ir++)
            {
                VMB[iq][idx[0]] = DENSS[ip][ir][2];
                iq += 1;
            }
        }
    }

    for (iq=0; iq<N_ACC; iq++)
    {
        if (fabs(VMB[iq][idx[0]]) > 10.)
        {
            if(debug==1)fprintf(myfile," iq VMB %i %e\n", iq, VMB[iq][idx[0]]);
            if(debug==1)fprintf(myfile,"ERROR");
            // exit(0);
            return 999;

        }
        FMB[iq][idx[DIM_LA]] = VMB[iq][idx[0]]-VMB[iq][idx[DIM_LA]];
    }

    if (METHOD == 0)   /* simple mixing */
    {
        relax = 0.4;
        for (iq=0; iq<N_ACC; iq++)
        {
            tmpm = relax*FMB[iq][idx[DIM_LA]];
            VMB[iq][idx[0]] = VMB[iq][idx[DIM_LA]]+tmpm;
            /*VMB[iq][idx[0]] = VMB[iq][idx[DIM_LA]]+relax*FMB[iq][idx[DIM_LA]];*/
        }
    }
    else   /* modified Broyden's method */
    {
        /* A. Baran et al., Phys. Rev. C 78 (2008) 014318 */
        relax = 0.6;
        if (M_ACC > 1)
        {
            itmp = DIM_LA-1;
            tmpf = tmpm = 0.;
            for (iq=0; iq<N_ACC; iq++)
            {
                DFMB[iq][idx[itmp]] = FMB[iq][idx[itmp+1]]-FMB[iq][idx[itmp]];
                DVMB[iq][idx[itmp]] = VMB[iq][idx[itmp+1]]-VMB[iq][idx[itmp]];
                tmpm += fabs(VMB[iq][idx[0]]);
                tmpf += QUAD(DFMB[iq][idx[itmp]]);
            }
            /*printf(" tmpm tmpf %e %e\n", tmpm, tmpf);*/
            if (tmpf > 0.) tmpf = 1./sqrt(tmpf);
            for (iq=0; iq<N_ACC; iq++)
            {
                DFMB[iq][idx[itmp]] *= tmpf;
                DVMB[iq][idx[itmp]] *= tmpf;
                UMB[iq][idx[itmp]] = relax*DFMB[iq][idx[itmp]]+DVMB[iq][idx[itmp]];
            }

            for (ik=0; ik<DIM_LA; ik++)
            {
                for (in=0; in<DIM_LA; in++)
                {
                    MAT[ik][in] = 0.;
                    for (iq=0; iq<N_ACC; iq++)
                        MAT[ik][in] += DFMB[iq][idx[in]]*DFMB[iq][idx[ik]];
                }
            }
            for (ik=0; ik<DIM_LA; ik++) MAT[ik][ik] += 0.0001;

            get_mat_inv(DIM_LA);

            for (ik=0; ik<DIM_LA; ik++)
            {
                CMB[ik] = 0.;
                for (iq=0; iq<N_ACC; iq++)
                    CMB[ik] += DFMB[iq][idx[ik]]*FMB[iq][idx[DIM_LA]];
            }

            for (in=0; in<DIM_LA; in++)
            {
                GMB[in] = 0.;
                for (ik=0; ik<DIM_LA; ik++)
                    GMB[in] += CMB[ik]*MAT[ik][in];
            }

            for (iq=0; iq<N_ACC; iq++)
            {
                for (ip=0; ip<DIM_LA; ip++)
                    GUMB[iq][idx[ip]] = UMB[iq][idx[ip]]*GMB[ip];
            }

        }

        for (iq=0; iq<N_ACC; iq++)
        {
            VMB[iq][idx[0]] = VMB[iq][idx[DIM_LA]]+relax*FMB[iq][idx[DIM_LA]];
        }

        for (iq=0; iq<N_ACC; iq++)
        {
            for (ip=0; ip<DIM_LA; ip++)
                VMB[iq][idx[0]] -= GUMB[iq][idx[ip]];
        }

    }

    /*
    iq = 0;
    for (ip=0; ip<N_PART; ip++) {
      if (IN_PART[ip] == 1) {
        for (ir=0; ir<NR; ir++) {
    DENS[ip][ir][2] = VMB[iq][idx[0]];
    iq += 1;
        }
        for (ir=0; ir<NR; ir++) {
    DENSS[ip][ir][2] = VMB[iq][idx[0]];
    iq += 1;
        }
      }
    }
    */

    iq = 0;
    dmax = 0.;
    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            for (ir=0; ir<NR; ir++)
            {
                tmpm = fabs(DENS[ip][ir][2] - VMB[iq][idx[0]]);
                if (tmpm > dmax) dmax = tmpm;
                DENS[ip][ir][2] = VMB[iq][idx[0]];
                iq += 1;
            }
            for (ir=0; ir<NR; ir++)
            {
                tmpm = fabs(DENS[ip][ir][2] - VMB[iq][idx[0]]);
                if (tmpm > dmax) dmax = tmpm;
                DENSS[ip][ir][2] = VMB[iq][idx[0]];
                iq += 1;
            }
        }
    }
    /*printf(" dmax %e\n", dmax);
      if (dmax > 1.) exit(0);*/

    return 0;
}
/*****************************************************************************/
int get_acc0(void)
{
    int ip,iq,ir,idx[REC_ACC];
    double relax,tmpm,dmax;

    for (ip=0; ip<REC_ACC; ip++)
    {
        idx[ip] = (ip+M_ACC0)%REC_ACC;
    }

    /*
    for (ip=0; ip<REC_ACC; ip++) {
      if ((idx[ip] < 0) || (idx[ip] > (REC_ACC+1))) {
       if(debug==1)fprintf(myfile," ip idx %i %i\n", idx[ip]);
        call_error(3333);
      }
    }
    */

    M_ACC0 += 1;

    iq = 0;
    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {

            if ((IPHASE > 0) && (PARTICLE[ip].rms > 0.))
            {
                /* vector density */
                for (ir=0; ir<NR; ir++)
                {
                    DENS_FOLD[ir] = DENS[ip][ir][2];
                }
                fold_dens(ip);
                for (ir=0; ir<NR; ir++)
                {
                    DENS[ip][ir][2] = DENS_FOLD[ir];
                }
                /* scalar density */
                for (ir=0; ir<NR; ir++)
                {
                    DENS_FOLD[ir] = DENSS[ip][ir][2];
                }
                fold_dens(ip);
                for (ir=0; ir<NR; ir++)
                {
                    DENSS[ip][ir][2] = DENS_FOLD[ir];
                }
            }

            for (ir=0; ir<NR; ir++)
            {
                VMB0[iq][idx[0]] = DENS[ip][ir][2];
                iq += 1;
            }
            for (ir=0; ir<NR; ir++)
            {
                VMB0[iq][idx[0]] = DENSS[ip][ir][2];
                iq += 1;
            }
        }
    }

    for (iq=0; iq<N_ACC0; iq++)
    {
        if (fabs(VMB[iq][idx[0]]) > 10.)
        {
            if(debug==1)fprintf(myfile," iq VMB %i %e\n", iq, VMB0[iq][idx[0]]);
            //  exit(0);
            return 999;
        }
        FMB0[iq][idx[DIM_LA]] = VMB0[iq][idx[0]]-VMB0[iq][idx[DIM_LA]];
    }

    /*
    for (iq=0; iq<N_ACC0; iq++) {
     if(debug==1)fprintf(myfile," %e", FMB0[iq][idx[DIM_LA]]);
    }
    if(debug==1)fprintf(myfile,"\n");
    */

    relax = 0.25;
    for (iq=0; iq<N_ACC0; iq++)
    {
        tmpm = relax*FMB0[iq][idx[DIM_LA]];
        VMB0[iq][idx[0]] = VMB0[iq][idx[DIM_LA]]+tmpm;
        /*VMB[iq][idx[0]] = VMB[iq][idx[DIM_LA]]+relax*FMB[iq][idx[DIM_LA]];*/
    }

    iq = 0;
    dmax = 0.;
    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            for (ir=0; ir<NR; ir++)
            {
                tmpm = fabs(DENS[ip][ir][2] - VMB0[iq][idx[0]]);
                if (tmpm > dmax) dmax = tmpm;
                DENS[ip][ir][2] = VMB0[iq][idx[0]];
                iq += 1;
            }
            for (ir=0; ir<NR; ir++)
            {
                tmpm = fabs(DENS[ip][ir][2] - VMB0[iq][idx[0]]);
                if (tmpm > dmax) dmax = tmpm;
                DENSS[ip][ir][2] = VMB0[iq][idx[0]];
                iq += 1;
            }
        }
    }
    /*printf(" dmax %e\n", dmax);
      if (dmax > 1.) exit(0);*/

    return 0;
}
/*****************************************************************************/
int init_acc(void)
{
    int ip,ir,iq;

    if (REC_ACC != (DIM_LA+1)) call_error(300);

    M_ACC = 0;
    for (ip=0; ip<DIM_LA; ip++)
    {
        for (iq=0; iq<DIM_LA; iq++)
        {
            AMB[ip][iq] = 0.;
        }
    }

    iq = 0;
    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            for (ir=0; ir<NR; ir++)
            {
                VMB[iq][DIM_LA] = DENS[ip][ir][2];
                /*printf(" ip iq %i %i\n", ip, iq);*/
                iq += 1;
            }
            for (ir=0; ir<NR; ir++)
            {
                VMB[iq][DIM_LA] = DENSS[ip][ir][2];
                /*printf(" ip iq %i %i\n", ip, iq);*/
                iq += 1;
            }
        }
    }
    N_ACC = iq;
    /*printf(" N_ACC = %i\n", N_ACC);*/

    for (iq=0; iq<N_ACC; iq++)
    {
        for (ip=0; ip<REC_ACC; ip++)
        {
            FMB[iq][ip] = DFMB[iq][ip] = DVMB[iq][ip] =
                                             UMB[iq][ip] = GUMB[iq][ip] = 0.;
        }
    }

    return 0;
}
/*****************************************************************************/
int init_acc0(void)
{
    int ip,ir,iq;

    if (REC_ACC != (DIM_LA+1)) call_error(300);

    M_ACC0 = 0;

    iq = 0;
    for (ip=0; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            for (ir=0; ir<NR; ir++)
            {
                VMB0[iq][DIM_LA] = DENS[ip][ir][2];
                /*printf(" ip iq %i %i\n", ip, iq);*/
                iq += 1;
            }
            for (ir=0; ir<NR; ir++)
            {
                VMB0[iq][DIM_LA] = DENSS[ip][ir][2];
                /*printf(" ip iq %i %i\n", ip, iq);*/
                iq += 1;
            }
        }
    }
    N_ACC0 = iq;
    /*printf(" N_ACC = %i\n", N_ACC);*/

    for (iq=0; iq<N_ACC0; iq++)
    {
        for (ip=0; ip<REC_ACC; ip++)
        {
            FMB0[iq][ip] = 0.;
        }
    }

    return 0;
}
/*****************************************************************************/
int get_mat_inv(int dim)
{
    int i1,i2;

    for (i1=0; i1<dim; i1++)
    {
        for (i2=0; i2<dim; i2++)
        {
            MAT2[i1][i2] = MAT[i1][i2];
        }
    }

    lu_decomposition(dim);

    for (i1=0; i1<dim; i1++)
    {
        for (i2=0; i2<dim; i2++) YVEC[i2] = YVEC2[i2] = 0.;
        YVEC[i1] = YVEC2[i1] = 1.;
        lu_backsubstitution(dim);
        improve(dim);
        for (i2=0; i2<dim; i2++) MAT3[i1][i2] = XVEC[i2];
    }

    for (i1=0; i1<dim; i1++)
    {
        for (i2=0; i2<dim; i2++)
        {
            MAT[i1][i2] = MAT3[i1][i2];
        }
    }


    return 0;
}
/*****************************************************************************/
int improve(int n)
{
    /* improves solution x of linear equations AX=Y
       adapted from page 48 in
       W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
       Numerical Recipes, Cambridge University Press, 1992 */
    int i, j;
    double sdp;

    for (i=0; i<n; i++)
    {
        XVEC[i] = YVEC[i];
    }

    for (i=0; i<n; i++)
    {
        sdp = -YVEC2[i];
        for (j=0; j<n; j++)
        {
            sdp += MAT2[i][j]*XVEC[j];
        }
        YVEC[i] = sdp;
    }

    lu_backsubstitution(n);

    for (i=0; i<n; i++)
    {
        XVEC[i] -= YVEC[i];
    }

    return 0;
}

/*****************************************************************************/
int lu_backsubstitution(int n)
{
    /* solves n linear equations AX=Y with LU decomposed NxN matrix A
       adapted from page 39 in
       W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
       Numerical Recipes, Cambridge University Press, 1992*/
    int i, i2, ii, j, ll;
    double sum;

    ii = -1;

    for (i=0; i<n; i++)
    {
        ll = INDX[i];
        sum = YVEC[ll];
        YVEC[ll] = YVEC[i];
        if (ii != -1)
        {
            for (j=ii; j<i; j++)
            {
                sum -= MAT[i][j]*YVEC[j];
            }
        }
        else
        {
            if (sum != 0.)
            {
                ii = i;
            }
        }
        YVEC[i] = sum;
    }

    for (i2=0; i2<n; i2++)
    {
        i = n-1-i2;
        sum = YVEC[i];
        for (j=i+1; j<n; j++)
        {
            sum -= MAT[i][j]*YVEC[j];
        }
        YVEC[i] = sum/MAT[i][i];
    }

    return 0;
}
/*****************************************************************************/
int lu_decomposition(int n)
{
    /* LU decomposition of NxN matrix A
       adapted from pages 38-39 in
       W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
       Numerical Recipes, Cambridge University Press, 1992 */
    int i, j, k, imax;
    double d, tiny, aamax, sum, dum, tmp;

    d = 1.;
    tiny = 1.e-20;

    for (i=0; i<n; i++)
    {
        aamax = 0.;
        for (j=0; j<n; j++)
        {
            tmp = fabs(MAT[i][j]);
            if (tmp > aamax) aamax=tmp;
        }
        if (aamax == 0.)
        {
            if(debug==1)fprintf(myfile," n M_ACC %i %i\n", n, M_ACC);
            for (j=0; j<4; j++)if(debug==1)fprintf(myfile," j idx[j] %i %i\n", j, INDEX[j]);
            if(debug==1)fprintf(myfile," MAT\n");
            for (j=0; j<3; j++)
            {
                for (k=0; k<3; k++)if(debug==1)fprintf(myfile," %e", MAT[j][k]);
                if(debug==1)fprintf(myfile,"\n");
            }
            if(debug==1)fprintf(myfile," FMB\n");
            for (j=0; j<N_ACC; j++)
            {
                for (k=0; k<REC_ACC; k++)if(debug==1)fprintf(myfile," %e", DFMB[j][INDEX[k]]);
                if(debug==1)fprintf(myfile,"\n");
            }
            /*
            for (j=0; j<N_CL; j++) {
            if(debug==1)fprintf(myfile," j = %i\n", j);
            if(debug==1)fprintf(myfile," BE_CL  %e\n", BE_CL[j]*HBARC);
            if(debug==1)fprintf(myfile," DBE_CL %e\n", DBE_CL[0][j]*HBARC);
            if(debug==1)fprintf(myfile," DBE_CL %e\n", DBE_CL[1][j]*HBARC);
            if(debug==1)fprintf(myfile," DBE_CL %e\n", DBE_CL[2][j]*HBARC);
            if(debug==1)fprintf(myfile," DBE_CL %e\n", DBE_CL[3][j]*HBARC);
            if(debug==1)fprintf(myfile," DEJ DEJT %e %e\n", DEJ, DEJT);
            }
            */
            call_error(3);
        }
        LUVV[i] = 1./aamax;
    }

    for (j=0; j<n; j++)
    {
        for (i=0; i<j; i++)
        {
            sum = MAT[i][j];
            for (k=0; k<i; k++)
            {
                sum -= MAT[i][k]*MAT[k][j];
            }
            MAT[i][j] = sum;
        }
        aamax = 0.;
        for (i=j; i<n; i++)
        {
            sum = MAT[i][j];
            for (k=0; k<j; k++)
            {
                sum -= MAT[i][k]*MAT[k][j];
            }
            MAT[i][j] = sum;
            dum = LUVV[i]*fabs(sum);
            if (dum >= aamax)
            {
                imax = i;
                aamax = dum;
            }
        }
        if (j != imax)
        {
            for (k=0; k<n; k++)
            {
                dum = MAT[imax][k];
                MAT[imax][k] = MAT[j][k];
                MAT[j][k] = dum;
            }
            d *= -1.;
            LUVV[imax] = LUVV[j];
        }
        INDX[j] = imax;
        if (MAT[j][j] == 0.) MAT[j][j] = tiny;
        if (j != n)
        {
            dum = 1./MAT[j][j];
            for (i=j+1; i<n; i++)
            {
                MAT[i][j] *= dum;
            }
        }
    }

    return 0;
}
/*****************************************************************************/
int get_chemical_potential_bar(double mu_b,double mu_q)
{
    double dmu_q,dmu_b;
    int ip;

    /*
    if(debug==1)fprintf(myfile," mu_b mu_q %e %e\n", mu_b*HBARC, mu_q*HBARC);
    get_nb_yq(mu_b,mu_q,1);
    for (ip=0; ip<N_PART; ip++) {
      if (IN_PART[ip] == 1) {
       if(debug==1)fprintf(myfile," ip S V rho %i %e %e %e\n",
         ip, SS[ip][0]*HBARC, VV[ip][0]*HBARC, DENS[ip][0][2]);
      }
      }*/


    if (ICOUNT == 3)
    {
        if(debug==1)fprintf(myfile," mu_b mu_q %e %e\n", mu_b*HBARC, mu_q*HBARC);
        dmu_b = 0.1;
        dmu_q = 0.05;
        for (ip=0; ip<1; ip++)
        {
            get_nb_yq(mu_b+dmu_b*(double)ip,mu_q,1);
            if(debug==1)fprintf(myfile," %e %e %e %e %i\n",
                                    (mu_b+dmu_b*(double)ip)*HBARC, mu_q*HBARC, F_AA, F_YQ, XCL);
        }
        //exit(0);
        return 0;
    }

    CI_NBYQ = 0;

    dmu_b = 0.1;
    dmu_q = 0.05;
    if(fit_asymmetry(mu_b,dmu_b,mu_q,dmu_q) == 999) return 999;

    mu_b = F_MU_B;
    mu_q = F_MU_Q;
    get_nb_yq(mu_b,mu_q,1);

    /*
    if(debug==1)fprintf(myfile," mu_b mu_q %e %e\n", mu_b*HBARC, mu_q*HBARC);
    for (ip=0; ip<N_PART; ip++) {
      if (IN_PART[ip] == 1) {
       if(debug==1)fprintf(myfile," ip S V rho %i %e %e %e\n",
         ip, SS[ip][0]*HBARC, VV[ip][0]*HBARC, DENS[ip][0][2]);
      }
    }
    for (ip=0; ip<N_AME11; ip++) {
      if (DENS_TON[ip] > 0.) {
        if(debug==1)fprintf(FILE_PLOT," %i %i %i %e\n",
          ip, NUC[ip].a, NUC[ip].z, DENS_TON[ip]);
      }
    }
    exit(0);
    */

    return 0;
}
/*****************************************************************************/
int fit_asymmetry(double mu_b,double dmu_b,double mu_q,double dmu_q)
{
    double xm,x0,xp,fm,f0,fp,dx,x0old,a,b,c,tmp1,tmp2,x1,x2,tmp3;
    int iwr,ic,ic0;

    iwr = 0;
    ic = ic0 = 0;

    /*if (ICOUNT == 19) iwr = 1;*/

    if (0 == 1)
    {
        for (ic=-20; ic<21; ic++)
        {
            x0 = (-1.3+0.01*(double)ic)/HBARC;
            fit_density(mu_b,dmu_b,x0);
            f0 = F_YQ;
            mu_b = F_MU_B;
            if(debug==1)fprintf(FILE_PLOT," %e %e %e\n", x0*HBARC, f0, mu_q*HBARC);
        }
        //exit(0);
        return 0;
    }

    dx = 2.*dmu_q;
    /*dx /= 100.;*/

    x0 = mu_q-0.5*dx;
    ic0 += 1;
    fit_density(mu_b,dmu_b,x0);
    f0 = F_YQ;
    mu_b = F_MU_B;
    if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);

    if (f0 > 0.)
    {
        do
        {
            xp = x0;
            fp = f0;
            x0 -= dx;
            ic0 += 1;
            fit_density(mu_b,dmu_b,x0);
            f0 = F_YQ;
            mu_b = F_MU_B;
            if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);
            /*if ((ICOUNT == 19) && (ic0 > 30)) exit(0);*/
            if(ic0 > 100) return 999; /* Stuck in an infinite loop here for some cases depending on the compiler.
                                      * Need to check thoroughly.*/
        }
        while (f0 > 0.);
        xm = x0;
        fm = f0;
    }
    else
    {
        if (f0 < 0.)
        {
            do
            {
                xm = x0;
                fm = f0;
                x0 += dx;
                ic0 += 1;
                fit_density(mu_b,dmu_b,x0);
                f0 = F_YQ;
                mu_b = F_MU_B;
                if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);
                if(ic0 > 100) return 999;
            }
            while (f0 < 0.);
            xp = x0;
            fp = f0;
        }
        else
        {
            F_MU_Q = x0;
            return 0;
        }
    }

    if (iwr == 1)if(debug==1)fprintf(myfile," xm fm %f %e\n", xm*HBARC, fm);
    if (iwr == 1)if(debug==1)fprintf(myfile," xp fp %f %e\n", xp*HBARC, fp);

    /*if (ICOUNT == 19) exit(0);*/

    /*
      x0 = 0.5*(xm+xp);
      fit_density(mu_b,dmu_b,x0);
      f0 = F_YQ;
      ic += 1;
      mu_b = F_MU_B;
    */

    x0old = xp;
    do
    {

        /*if (iwr == 1)if(debug==1)fprintf(myfile,"\n xm fm %f %e\n", xm*HBARC, fm);
          if (iwr == 1)if(debug==1)fprintf(myfile," xp fp %f %e\n\n", xp*HBARC, fp);*/

        x0 = 0.5*(xm+xp);
        /*x0 = (fp*xm-fm*xp)/(fp-fm);*/
        fit_density(mu_b,dmu_b,x0);
        f0 = F_YQ;
        ic += 1;
        mu_b = F_MU_B;
        if (iwr == 1)if(debug==1)fprintf(myfile," ic x0 f0 %i %f %e\n", ic, x0*HBARC, f0);
        a = f0;
        b = 0.5*(fp-fm);
        c = 0.5*(fp-2.*f0+fm);
        tmp1 = QUAD(b);
        tmp2 = 4.*a*c;
        tmp3 = 0.5*(xp-xm);
        if (tmp1 > tmp2)
        {
            if (fabs(tmp2/tmp1) > 1.e-06)
            {
                tmp1 -= tmp2;
                tmp1 = sqrt(tmp1);
                x1 = (-b+tmp1)/(2.*c);
                x2 = (-b-tmp1)/(2.*c);
            }
            else
            {
                x1 = x2 = -a/b;
            }
            if (fabs(x2) < 1.) x1 = x2;
            if (fabs(x1) < 1.)
            {
                if (fabs(x1) < 0.25)
                {
                    x2 = x0+2.*x1*tmp3;
                }
                else
                {
                    x2 = x0;
                }
                if (fabs(x1) > 0.75)
                {
                    if (x1 > 0.) x2 = x0+(2.*x1-1.)*tmp3;
                    if (x1 < 0.) x2 = x0+(2.*x1+1.)*tmp3;
                }
                x1 = x0+x1*tmp3;

                if (f0 < 0.)
                {
                    xm = x0;
                    fm = f0;
                }
                else
                {
                    xp = x0;
                    fp = f0;
                }
                if ((xm < x1) && (x1 < xp))
                {
                    x0 = x1;
                    fit_density(mu_b,dmu_b,x0);
                    f0 = F_YQ;
                    ic += 1;
                    mu_b = F_MU_B;
                    if (iwr == 1)if(debug==1)fprintf(myfile," ic x1 f1 %i %f %e\n", ic, x0*HBARC, f0);
                    if (f0 < 0.)
                    {
                        xm = x0;
                        fm = f0;
                    }
                    else
                    {
                        xp = x0;
                        fp = f0;
                    }
                }
                if ((xm < x2) && (x2 < xp))
                {
                    x0 = x2;
                    fit_density(mu_b,dmu_b,x0);
                    f0 = F_YQ;
                    ic += 1;
                    mu_b = F_MU_B;
                    if (iwr == 1)if(debug==1)fprintf(myfile," ic x2 f2 %i %f %e\n", ic, x0*HBARC, f0);
                    if (f0 < 0.)
                    {
                        xm = x0;
                        fm = f0;
                    }
                    else
                    {
                        xp = x0;
                        fp = f0;
                    }
                }
            }
        }


        if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);

        if ((xm < x0) && (x0 < xp))
        {
            if (f0 < 0.)
            {
                xm = x0;
                fm = f0;
            }
            else
            {
                xp = x0;
                fp = f0;
            }
        }

        dx = fabs(x0-x0old);
        x0old = x0;

        /*dx = xp-xm;*/
        if (iwr == 1)if(debug==1)fprintf(myfile," dx xm xp fm fp %e %f %f %e %e\n",
                                                 dx, xm*HBARC, xp*HBARC, fm, fp);

    }
    while ((dx > 0.5e-09) && (ic < 200));
    /*} while (((fabs(f0) > 1.e-08) || (dx > 0.5e-09)) && (ic < 200));*/

    CI_YQ = ic; /*+ic0;*/

    F_MU_Q = x0;

    return 0;
}
/*****************************************************************************/
int fit_density(double mu_b,double dmu_b,double mu_q)
{
    double xm,x0,xp,fm,f0,fp,dx,a,b,c,tmp1,tmp2,x1,x2,tmp3,x0old,
           mu_b_max,mu_b_min;
    int iwr,ic,ic0,icon;

    iwr = 0;
    ic = ic0 = icon = 0;
    /*if (ICOUNT == 19) iwr = 1;*/

    /* test */
    /*
    if (XCL == 1) {
    for (ic=0; ic<10; ic++) {
      x0 = (920.+1.*(double)ic)/HBARC;
      get_nb_yq(x0,mu_q,0);
      f0 = F_AA;
     if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);
    }
    exit(0);
    }
    */

    mu_b_min = -1.e+99;
    mu_b_max =  1.e+99;

    /* limits for condensation */

    if (mu_b > mu_b_max) mu_b = mu_b_max;
    if (mu_b < mu_b_min) mu_b = mu_b_min;

    dx = 2.*dmu_b;
    if (dx > (mu_b_max-mu_b_min)) dx = 0.25*(mu_b_max-mu_b_min);

    if (IN_CL == 1) dx = 1./HBARC;

    x0 = mu_b;
    ic0 += 1;
    get_nb_yq(x0,mu_q,0);
    f0 = F_AA;
    if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);

    if (f0 >= 0.)
    {
        do
        {
            xp = x0;
            fp = f0;
            x0 -= dx;
            if (x0 < mu_b_min) x0 = mu_b_min;
            ic0 += 1;
            get_nb_yq(x0,mu_q,0);
            f0 = F_AA;
            if (iwr == 1)if(debug==1)fprintf(myfile,"+ x0 f0 %f %e\n", x0*HBARC, f0);
            if(ic0 > 200) return 999; /* Stuck in an infinite loop here for some cases depending on the compiler.
                                      * Need to check thoroughly.*/
        }
        while (f0 > 0.);
        xm = x0;
        fm = f0;
    }
    else
    {
        if (f0 < 0.)
        {
            do
            {
                xm = x0;
                fm = f0;
                x0 += dx;
                if (x0 > mu_b_max) x0 = mu_b_max;
                ic0 += 1;
                get_nb_yq(x0,mu_q,0);
                f0 = F_AA;
                if (iwr == 1)if(debug==1)fprintf(myfile,"- x0 f0 %f %e\n", x0*HBARC, f0);
                if(ic0 > 200) return 999; /* Stuck in an infinite loop here for some cases depending on the compiler.
                                      * Need to check thoroughly.*/
            }
            while (f0 < 0.);
            xp = x0;
            fp = f0;
        }
        else
        {
            F_MU_B = x0;
            return 0;
        }
    }


    x0old = xp;
    do
    {
        x0 = 0.5*(xm+xp);
        /*x0 = (fp*xm-fm*xp)/(fp-fm);*/
        get_nb_yq(x0,mu_q,0);
        ic += 1;
        f0 = F_AA;
        if (iwr == 1)if(debug==1)fprintf(myfile," ic x0 f0 %i %f %e\n", ic, x0*HBARC, f0);
        a = f0;
        b = 0.5*(fp-fm);
        c = 0.5*(fp-2.*f0+fm);
        tmp1 = QUAD(b);
        tmp2 = 4.*a*c;
        tmp3 = 0.5*(xp-xm);
        if (tmp1 > tmp2)
        {
            if (fabs(tmp2/tmp1) > 1.e-06)
            {
                tmp1 -= tmp2;
                tmp1 = sqrt(tmp1);
                x1 = (-b+tmp1)/(2.*c);
                x2 = (-b-tmp1)/(2.*c);
            }
            else
            {
                x1 = x2 = -a/b;
            }
            if (fabs(x2) < 1.) x1 = x2;
            if (fabs(x1) < 1.)
            {
                if (fabs(x1) < 0.25)
                {
                    x2 = x0+2.*x1*tmp3;
                }
                else
                {
                    x2 = x0;
                }
                if (fabs(x1) > 0.75)
                {
                    if (x1 > 0.) x2 = x0+(2.*x1-1.)*tmp3;
                    if (x1 < 0.) x2 = x0+(2.*x1+1.)*tmp3;
                }
                x1 = x0+x1*tmp3;

                if (f0 < 0.)
                {
                    xm = x0;
                    fm = f0;
                }
                else
                {
                    xp = x0;
                    fp = f0;
                }
                if ((xm < x1) && (x1 < xp))
                {
                    x0 = x1;
                    get_nb_yq(x0,mu_q,0);
                    ic += 1;
                    f0 = F_AA;
                    if (iwr == 1)if(debug==1)fprintf(myfile," ic x1 f1 %i %f %e\n", ic, x0*HBARC, f0);
                    if (f0 < 0.)
                    {
                        xm = x0;
                        fm = f0;
                    }
                    else
                    {
                        xp = x0;
                        fp = f0;
                    }
                }
                if ((xm < x2) && (x2 < xp))
                {
                    x0 = x2;
                    get_nb_yq(x0,mu_q,0);
                    ic += 1;
                    f0 = F_AA;
                    if (iwr == 1)if(debug==1)fprintf(myfile," ic x2 f2 %i %f %e\n", ic, x0*HBARC, f0);
                    if (f0 < 0.)
                    {
                        xm = x0;
                        fm = f0;
                    }
                    else
                    {
                        xp = x0;
                        fp = f0;
                    }
                }
            }
        }


        if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);

        if ((xm < x0) && (x0 < xp))
        {
            if (f0 < 0.)
            {
                xm = x0;
                fm = f0;
            }
            else
            {
                xp = x0;
                fp = f0;
            }
        }

        dx = fabs(x0-x0old);
        x0old = x0;

        /*dx = xp-xm;*/
        if (iwr == 1)if(debug==1)fprintf(myfile," dx xm xp fm fp %e %f %f %e %e\n",
                                                 dx, xm*HBARC, xp*HBARC, fm, fp);

        /*if (ic > 100) exit(0);*/

        /*} while ((dx > 0.5e-09) && (ic < 200));*/
    }
    while (((fabs(f0) > 1.e-08) || (dx > 0.5e-09)) && (ic < 200));


    CI_NB = ic;
    /*+ic0;*/

    F_MU_B = x0;

    return 0;
}
/*****************************************************************************/
int get_nb_yq(double mu_b,double mu_q,int ic)
{
    double mu,tmpb,tmpq,rho;
    int ip,i1;

    /* proton and neutron chemical potentials*/
    PARTICLE[0].mu = mu_b+mu_q;
    PARTICLE[1].mu = mu_b;

    if (ICOUNT == 3)
    {
        if(debug==1)fprintf(myfile," mu_b %e\n", mu_b*HBARC);
        if(debug==1)fprintf(myfile," mu_q %e\n", mu_q*HBARC);
        if(debug==1)fprintf(myfile," p mu %e\n", PARTICLE[0].mu*HBARC);
        if(debug==1)fprintf(myfile," n mu %e\n", PARTICLE[1].mu*HBARC);
    }

    /* NSE with nonrelativistic Maxwell-Boltzmann distribution ? */

    /* RMF */

    /* get_rho_part(double mu, int ip, int ic)*/

    tmpb = tmpq = 0.;

    /* nucleons */
    for (ip=0; ip<2; ip++)
    {
        mu = PARTICLE[ip].mu;
        rho = get_rho_part(mu,ip,ic);
        /*printf(" ip mu rho %i %f %e\n", ip, mu*HBARC, rho);*/
        tmpb += rho*(double)PARTICLE[ip].a;
        tmpq += rho*(double)PARTICLE[ip].z;
    }
    /*if (ICOUNT == 5) tmpb = tmpq = 0.;*/

    /* other particles */
    /* what about bosons = antibosons? */

    /* clusters */
    if (XCL == 1)
    {
        for (i1=0; i1<N_CL; i1++)
        {
            ip = 9+i1;
            if (IN_PART[ip] == 1)
            {
                mu = PARTICLE[ip].mu = PARTICLE[0].mu*(double)PARTICLE[ip].z
                                       +PARTICLE[1].mu*(double)PARTICLE[ip].n;
                rho = get_rho_part(mu,ip,ic);
                /*printf("d %f %e\n", mu*HBARC,rho);*/
                /*printf(" mu S rho %e %e %e\n", mu*HBARC, SS[ip][0]*HBARC, rho);*/
                if (ICOUNT == 3)
                {
                    /*printf(" ip rho m s v mu %i %e %e %e %e %e\n",
                      ip, rho, PARTICLE[ip].m*HBARC,
                      SS[ip][0]*HBARC, VV[ip][0]*HBARC, mu*HBARC);*/
                    if(debug==1)fprintf(myfile," ip g %i %f\n", ip, PARTICLE[ip].g);
                    if(debug==1)fprintf(myfile," ip rho mu m-s mu-v %i %e %e %e %e\n",
                                            ip, rho, mu*HBARC,
                                            (PARTICLE[ip].m-SS[ip][0])*HBARC,
                                            (mu-VV[ip][0])*HBARC);
                    // exit(0);
                    return 0;

                }
                tmpb += rho*(double)PARTICLE[ip].a;
                tmpq += rho*(double)PARTICLE[ip].z;
            }
        }
    }

    /* hyperons */
    if (XHYP == 1)
    {
        if (IN_HYP == 1)
        {
            for (i1=0; i1<N_HYP; i1++)
            {
                ip = 9+N_CL+i1;
                if (IN_PART[ip] == 1)
                {
                    mu = PARTICLE[ip].mu = mu_q*(double)PARTICLE[ip].z
                                           +mu_b*(double)PARTICLE[ip].a;
                    rho = get_rho_part(mu,ip,ic);
                    tmpb += rho*(double)PARTICLE[ip].a;
                    tmpq += rho*(double)PARTICLE[ip].z;
                }
            }
        }
    }

    /* thermal mesons */
    if (IN_TMES == 1)
    {
        for (i1=0; i1<N_TMES; i1++)
        {
            ip = 9+N_CL+N_HYP+i1;
            if (IN_PART[ip] == 1)
            {
                mu = PARTICLE[ip].mu = mu_q*(double)PARTICLE[ip].z
                                       +mu_b*(double)PARTICLE[ip].a;
                rho = get_rho_part(mu,ip,ic);
                tmpb += rho*(double)PARTICLE[ip].a;
                tmpq += rho*(double)PARTICLE[ip].z;
            }
        }
    }

    /* table of nuclei */
    if (IN_TON == 1)
    {
        get_dens_ton(ic);
        tmpb += DENS_TON_B;
        tmpq += DENS_TON_Q;
    }

    /* condensation contribution */

    /* figures of merit */

    F_AA = (tmpb-AA_WS)/AA_WS;

    if (tmpb > 0.)
    {
        F_YQ = tmpq/tmpb-Y_Q;
    }
    else
    {
        F_YQ = tmpq/AA_WS-Y_Q;
    }

    CI_NBYQ += 1;

    return 0;
}
/*****************************************************************************/
int get_dens_ton(int ic)
{
    int ip,iq;
    double arg[N_AME11],fac,tmp;

    fac = pow(TT/TPI,1.5);

    /*printf("\n mu_p mu_n %f %f\n", PARTICLE[0].mu*HBARC, PARTICLE[1].mu*HBARC);*/

    for (ip=0; ip<N_NUC; ip++)
    {
        ISORT[ip] = ip;
        NUC[ip].ms = NUC[ip].m-SS_NUC[ip];
        NUC[ip].mu = PARTICLE[0].mu*(double)NUC[ip].z
                     +PARTICLE[1].mu*(double)NUC[ip].n;
        arg[ip] = (NUC[ip].mu-NUC[ip].ms-VV_NUC[ip])/TT;
        arg[ip] += 1.5*log(NUC[ip].ms);
        arg[ip] += log(NUC[ip].g);
        /*printf(" %i   %i %i   %f %f %f\n",
          ip, NUC[ip].a, NUC[ip].z, NUC[ip].mu*HBARC, NUC[ip].ms*HBARC, arg[ip]);*/
        VSORT[ip] = -arg[ip];
    }
    heapsort(N_NUC);

    DENS_TON_B = DENS_TON_Q = 0.;

    if (ic == 0)
    {
        for (ip=(N_NUC-1); ip>-1; ip--)
        {
            tmp = arg[ISORT[0]]-40.;
            iq = ISORT[ip];
            if (arg[iq] > tmp)
            {
                if (arg[iq] > 220.) arg[iq] = 220.;
                DENS_TON[iq] = fac*exp(arg[iq]);
                /*printf(" %i %i %e %e\n", ip, iq, arg[ISORT[iq]], DENS_TON[iq]);*/
                DENS_TON_B += DENS_TON[iq]*(double)NUC[iq].a;
                DENS_TON_Q += DENS_TON[iq]*(double)NUC[iq].z;
            }
            else
            {
                DENS_TON[iq] = 0.;
            }
        }
    }
    else
    {
        NP_NUC = 0;
        for (ip=(N_NUC-1); ip>-1; ip--)
        {
            tmp = arg[ISORT[0]]-40.;
            iq = ISORT[ip];
            if (arg[iq] > tmp)
            {
                NP_NUC += 1;
                DENS_TON[iq] = fac*exp(arg[iq]);
                DENSS_TON[iq] = DENS_TON[iq]*(1.-1.5*TT/NUC[iq].ms);
                /*printf(" %i %i %e %e %e\n",
                  ip, iq, arg[ISORT[iq]], DENS_TON[iq], DENSS_TON[iq]);*/
                DENS_TON_B += DENS_TON[iq]*(double)NUC[iq].a;
                DENS_TON_Q += DENS_TON[iq]*(double)NUC[iq].z;
            }
            else
            {
                DENS_TON[iq] = DENSS_TON[iq] = 0.;
            }
            /*printf(" dens_ton[0] %e\n", DENS_TON[0]);*/
        }
        /*exit(0);*/
    }
    /*
    if(debug==1)fprintf(myfile," %e %e\n", DENS_NUC_B, DENS_NUC_Q);
    exit(0);
    */

    return 0;
}
/*****************************************************************************/
int get_test()
{
    int ix;
    double mu,tmp,tmps;

    if(debug==1)fprintf(myfile,"hallo test\n");
    IPHASE = 0;
    VV[2][0] = 0.;
    SS[2][0] = 0.;
    TT = 75./HBARC;
    PARTICLE[2].st = 1;
    PARTICLE[2].m = 500./HBARC;
    for (ix=-100; ix<100; ix++)
    {
        mu = 6.*(double)ix/HBARC;
        PARTICLE[2].st = 1;
        tmp = get_rho_part(mu,2,1);
        PARTICLE[2].st = 11;
        tmps = get_rho_part(mu,2,1);
        /*tmps = DENSS[2][0][0]-DENSS[2][0][1];*/
        if(debug==1)fprintf(FILE_PLOT," %e %e %e %e\n",
                                mu*HBARC, tmp, tmps, fabs(tmps-tmp)/tmp);
    }

    return 0;
}
/*****************************************************************************/
int heapsort(int n)
{
    /* Sorting of vector by heapsort algorithm
       adapted from page 329 in
       W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery
       Numerical Recipes in Fortran, Cambridge University Press, 1992 */
    int i, ir, j, l, irra;
    double rra;

    if (n < 2) return 0;
    l = n/2+1;
    ir = n;
    while (ir > 0)
    {
        if (l > 1)
        {
            l -= 1;
            rra = VSORT[l-1];
            irra = ISORT[l-1];
        }
        else
        {
            ir -= 1;
            rra = VSORT[ir];
            irra = ISORT[ir];
            VSORT[ir] = VSORT[0];
            ISORT[ir] = ISORT[0];
            if (ir == 1)
            {
                VSORT[0] = rra;
                ISORT[0] = irra;
                return 0;
            }
        }
        i = l;
        j = l+l;
        while (j <= ir)
        {
            if (j < ir)
            {
                if (VSORT[j-1] < VSORT[j]) j += 1;
            }
            if (rra < VSORT[j-1])
            {
                VSORT[i-1] = VSORT[j-1];
                ISORT[i-1] = ISORT[j-1];
                i = j;
                j += j;
            }
            else
            {
                j = ir+1;
            }
        }
        VSORT[i-1] = rra;
        ISORT[i-1] = irra;
    }
    return 0;
}
/*****************************************************************************/
double get_rho_part(double mu,int ip,int ic)
{
    double rho,tmpv[DIM_R];
    int ir;

    F_ST = PARTICLE[ip].st;
    F_G  = PARTICLE[ip].g;
    F_M  = PARTICLE[ip].m;
    /*if (ip == 9) {
     if(debug==1)fprintf(myfile,"d m %f %f\n", PARTICLE[ip].m*HBARC, MM[ip][0]*HBARC);
      }*/

    if (ic == 0)
    {
        /* fitting */
        if (IPHASE > 0)
        {
            for (ir=0; ir<NR; ir++)
            {
                F_S = SS[ip][ir];
                F_V = VV[ip][ir];
                get_rho(mu);
                tmpv[ir]  = F_RHOP-F_RHOA;
            }
            rho = FPI*r_integrate(0,2,tmpv);
            if (IREXT == 1)
            {
                rho += VX*tmpv[NRP];
            }
        }
        else
        {
            F_S = SS[ip][0];
            F_V = VV[ip][0];
            get_rho(mu);
            rho = F_RHOP-F_RHOA;
        }
    }
    else
    {
        /* full calculation */
        if (IPHASE > 0)
        {
            for (ir=0; ir<NR; ir++)
            {
                F_S = SS[ip][ir];
                F_V = VV[ip][ir];
                get_rho(mu);
                DENS[ip][ir][0] = F_RHOP;
                DENS[ip][ir][1] = F_RHOA;
                DENS[ip][ir][2] = tmpv[ir] = F_RHOP-F_RHOA;
                get_rhos(mu);
                DENSS[ip][ir][0] = F_RHOP;
                DENSS[ip][ir][1] = F_RHOA;
                DENSS[ip][ir][2] = F_RHOP+F_RHOA;
            }
            rho = FPI*r_integrate(0,2,tmpv);
            if (IREXT == 1)
            {
                rho += VX*tmpv[NRP];
            }
        }
        else
        {
            F_S = SS[ip][0];
            F_V = VV[ip][0];
            get_rho(mu);
            DENS[ip][0][0] = F_RHOP;
            DENS[ip][0][1] = F_RHOA;
            DENS[ip][0][2] = rho = F_RHOP-F_RHOA;
            /*printf(" mu F_M F_S F_V TT rho %e %e %e %e %e %e %e\n",
               mu*HBARC, F_M*HBARC,  F_S*HBARC, F_V*HBARC, TT*HBARC, F_RHOP,
               PARTICLE[ip].g*exp((mu-F_M+F_S-F_V)/TT)/CUBE(sqrt(TPI/(TT*(F_M-F_S)))));*/
            get_rhos(mu);
            DENSS[ip][0][0] = F_RHOP;
            DENSS[ip][0][1] = F_RHOA;
            DENSS[ip][0][2] = F_RHOP+F_RHOA;
        }
    }
    /*
    if (rho > 100.) {
     if(debug==1)fprintf(myfile," ip rho m-S+V mu %i %e %e %e\n",
       ip, rho, (PARTICLE[ip].m-SS[ip][0]+VV[ip][0])*HBARC, mu*HBARC);
      exit(0);
    }
    */

    return rho;
}
/*****************************************************************************/
int get_rho(double mu)
{
    int ipa;
    double fac,ml,tmp,kf,res1,mup,tmp2;

    F_RHOP = F_RHOA = 0.;

    F_MD = F_M-F_S;
    F_MD2 = QUAD(F_MD);
    fac = F_G/TPI2;

    if (TT > 0.)
    {
        switch(F_ST)
        {
        case -1:   /* Bose-Einstein, particle = antiparticle */
        {
            /* limits on chemical potential, only mu = F_V possible? */
            mup = mu;
            if (mu > (F_V+fabs(F_MD))) mup = F_V+fabs(F_MD);
            IBX = ipa = 0;
            get_simint_bounds_b(mup,ipa);
            tmp = fac*get_simint(10,1.e-08,11,F_KM,F_KP,mup,ipa);
            if (ipa == 0)
            {
                F_RHOP = tmp;
            }
            else
            {
                F_RHOA = tmp;
            }
            break;
        }
        case 0:   /* Bose-Einstein, particle != antiparticle */
        {
            /*IBX = 0;
            if (fabs(F_V+F_MD-mu)/TT < 1.e-07) IBX =  1;
            if (fabs(F_V-F_MD-mu)/TT < 1.e-07) IBX = -1;*/
            /* limits on chemical potential */
            mup = mu;
            if (mu > (F_V+fabs(F_MD))) mup = F_V+fabs(F_MD);
            if (mu < (F_V-fabs(F_MD))) mup = F_V-fabs(F_MD);
            for (ipa=0; ipa<2; ipa++)
            {
                IBX = ipa;
                get_simint_bounds_b(mup,ipa);
                tmp  = fac*get_simint(10,1.e-08,11,F_KM,F_KP,mup,ipa);
                if (ipa == 0)
                {
                    F_RHOP = tmp;
                }
                else
                {
                    F_RHOA = tmp;
                }
            }
            break;
        }
        case 2:   /* Maxwell-Boltzmann, non-relativistic */
        {
            tmp = sqrt(fabs(F_MD)*TT/TPI);
            /*tmp = sqrt(fabs(F_M)*TT/TPI);*/
            tmp = CUBE(tmp);
            /*printf(" F_G tmp %e %e %e\n", F_G, tmp, exp((mu-fabs(F_MD)-F_V)/TT));*/
            F_RHOP = F_G*exp((mu-fabs(F_MD)-F_V)/TT)*tmp;
            F_RHOA = F_G*exp((F_V-fabs(F_MD)-mu)/TT)*tmp;
            /*printf(" P A %e %e\n", F_RHOP, F_RHOA);*/
            /* hallo */
            /*F_RHOA = 0.;*/
            break;
        }
        case 3:   /* Maxwell-Boltzmann, relativistic */
        {
            tmp = sqrt(fabs(F_MD)*TT/TPI);
            tmp = CUBE(tmp);
            F_RHOP = F_G*exp((mu-fabs(F_MD)-F_V)/TT)*tmp;
            F_RHOA = F_G*exp((F_V-fabs(F_MD)-mu)/TT)*tmp;
            tmp = get_k2_red(fabs(F_MD)/TT);
            F_RHOP *= tmp;
            F_RHOA *= tmp;
            break;
        }
        case 1:   /* Fermi-Dirac */
        {
            /* particle */
            ipa = 0;
            tmp2 = (mu-F_V-fabs(F_MD))/TT;
            if (tmp2 < -10.)
            {
                tmp = fabs(F_MD)/TT;
                tmp2 = exp(tmp2);
                F_RHOP = tmp2*get_k2_red(tmp);
                F_RHOP -= QUAD(tmp2)*get_k2_red(2.*tmp)/2.828427125;
                tmp = sqrt(fabs(F_MD)*TT/TPI);
                tmp = F_G*CUBE(tmp);
                F_RHOP *= tmp;
            }
            else
            {
                get_simint_bounds_f(mu,ipa);
                if (F_KP > F_KM)
                {
                    res1 = get_simint(5,1.e-08,1,F_KM,F_KP,mu,ipa);
                }
                else
                {
                    res1 = 0.;
                }
                res1 += CUBE(F_KM)/3.;
                F_RHOP = res1*fac;
            }
            /* antiparticle */
            ipa = 1;
            tmp2 = (F_V-mu-fabs(F_MD))/TT;
            if (tmp2 < -10.)
            {
                tmp = fabs(F_MD)/TT;
                tmp2 = exp(tmp2);
                F_RHOA = tmp2*get_k2_red(tmp);
                F_RHOA -= QUAD(tmp2)*get_k2_red(2.*tmp)/2.828427125;
                tmp = sqrt(fabs(F_MD)*TT/TPI);
                tmp = F_G*CUBE(tmp);
                F_RHOA *= tmp;
            }
            else
            {
                get_simint_bounds_f(mu,ipa);
                if (F_KP > F_KM)
                {
                    res1 = get_simint(5,1.e-08,1,F_KM,F_KP,mu,ipa);
                }
                else
                {
                    res1 = 0.;
                }
                res1 += CUBE(F_KM)/3.;
                F_RHOA = res1*fac;
            }
            break;
        }
        default:   /* Fermi-Dirac, original */
        {
            for (ipa=0; ipa<2; ipa++)
            {
                get_simint_bounds_f(mu,ipa);
                if (F_KP > F_KM)
                {
                    res1 = get_simint(5,1.e-08,1,F_KM,F_KP,mu,ipa);
                }
                else
                {
                    res1 = 0.;
                }
                res1 += CUBE(F_KM)/3.;
                if (ipa == 0)
                {
                    F_RHOP = res1;
                }
                else
                {
                    F_RHOA = res1;
                }
            }
            F_RHOP *= fac;
            F_RHOA *= fac;
        }
        }
    }
    else
    {
        if (F_ST == 1)   /* Fermi-Dirac */
        {
            ml = mu-F_V;
            tmp = QUAD(ml)-F_MD2;
            if (tmp > 0.)
            {
                kf = sqrt(tmp);
                tmp *= fac*kf/3.;
                if (ml >= 0.)
                {
                    F_RHOP = tmp;
                }
                else
                {
                    F_RHOA = tmp;
                }
            }
        }
    }

    /*printf(" F_RHOP F_RHOA %e %e\n", F_RHOP, F_RHOA);*/

    return 0;
}
/*****************************************************************************/
int get_rhos(double mu)
{
    int ipa;
    double fac,ml,tmp,kf,res1,mup;

    F_RHOP = F_RHOA = 0.;

    F_MD = F_M-F_S;
    F_MD2 = QUAD(F_MD);
    fac = F_G/TPI2;

    if (TT > 0.)
    {
        switch(F_ST)
        {
        case -1:   /* Bose-Einstein, particle = antiparticle */
        {
            /* limits on chemical potential, only mu = F_V = 0 possible? */
            mup = mu;
            if (mu > (F_V+fabs(F_MD))) mup = F_V+fabs(F_MD);
            for (ipa=0; ipa<1; ipa++)
            {
                IBX = ipa;
                get_simint_bounds_b(mup,ipa);
                tmp  = fac*get_simint(10,1.e-08,13,F_KM,F_KP,mup,ipa);
                if (ipa == 0)
                {
                    F_RHOP = tmp;
                }
                else
                {
                    F_RHOA = tmp;
                }
            }
            break;
            break ;
        }
        case 0:   /* Bose-Einstein, particle != antiparticle */
        {
            /*IBX = 0;
            if (fabs(F_V+F_MD-mu)/TT < 1.e-07) IBX =  1;
            if (fabs(F_V-F_MD-mu)/TT < 1.e-07) IBX = -1;*/
            /* limits on chemical potential */
            mup = mu;
            if (mu > (F_V+fabs(F_MD))) mup = F_V+fabs(F_MD);
            if (mu < (F_V-fabs(F_MD))) mup = F_V-fabs(F_MD);
            for (ipa=0; ipa<2; ipa++)
            {
                IBX = ipa;
                get_simint_bounds_b(mup,ipa);
                tmp  = fac*get_simint(10,1.e-08,13,F_KM,F_KP,mup,ipa);
                if (ipa == 0)
                {
                    F_RHOP = tmp;
                }
                else
                {
                    F_RHOA = tmp;
                }
            }
            break;
        }
        case 2:   /* Maxwell-Boltzmann, non-relativistic */
        {
            tmp = sqrt(fabs(F_MD)*TT/TPI);
            /*tmp = sqrt(fabs(F_M)*TT/TPI);*/
            tmp = CUBE(tmp)*(1.-1.5*TT/fabs(F_MD));
            F_RHOP = F_G*exp((mu-fabs(F_MD)-F_V)/TT)*tmp;
            F_RHOA = F_G*exp((F_V-fabs(F_MD)-mu)/TT)*tmp;
            break;
        }
        case 3:   /* Maxwell-Boltzmann, relativistic */
        {
            tmp = sqrt(fabs(F_MD)*TT/TPI);
            tmp = CUBE(tmp);
            F_RHOP = F_G*exp((mu-fabs(F_MD)-F_V)/TT)*tmp;
            F_RHOA = F_G*exp((F_V-fabs(F_MD)-mu)/TT)*tmp;
            tmp = get_k1_red(fabs(F_MD)/TT);
            F_RHOP *= tmp;
            F_RHOA *= tmp;
            break;
        }
        default:   /* Fermi-Dirac */
        {
            for (ipa=0; ipa<2; ipa++)
            {
                get_simint_bounds_f(mu,ipa);
                if (F_KP > F_KM)
                {
                    res1 = get_simint(5,1.e-08,3,F_KM,F_KP,mu,ipa);
                }
                else
                {
                    res1 = 0.;
                }
                if (F_KM > 0.)
                {
                    ml = sqrt(QUAD(F_KM)+F_MD2);
                    res1 += 0.5*fabs(F_MD)*(F_KM*ml-F_MD2*log((F_KM+ml)/fabs(F_MD)));
                }
                if (ipa == 0)
                {
                    F_RHOP = res1;
                }
                else
                {
                    F_RHOA = res1;
                }
            }
            F_RHOP *= fac;
            F_RHOA *= fac;
        }
        }
    }
    else
    {
        if (F_ST == 1)   /* Fermi-Dirac */
        {
            ml = mu-F_V;
            tmp = QUAD(ml)-F_MD2;
            if (tmp > 0.)
            {
                kf = sqrt(tmp);
                tmp = 0.5*fabs(F_MD)*
                      (kf*fabs(ml)-F_MD2*log((kf+fabs(ml))/fabs(F_MD)))*fac;
                if (ml >= 0.)
                {
                    F_RHOP = tmp;
                }
                else
                {
                    F_RHOA = tmp;
                }
            }
        }
    }

    return 0;
}
/*****************************************************************************/
double get_chemical_potential_lep(void)
{
    double fc,tmp,x0,f0,dx,xm,fm,xp,fp,x0old,a,b,c,tmp1,tmp2,tmp3,x1,x2;
    int iwr,ic0,ic;

    iwr = 0;

    ic0 = ic = 0;

    if (ZZ_WS != 0.)
    {
        fc = ZZ_WS;

        tmp = fc;
        if (IPHASE > 0) tmp /= V_WS;
        if (iwr == 1)if(debug==1)fprintf(myfile," ZZ_WS tmp %e %e\n", ZZ_WS, tmp);
        x0 = pow(3.*PI2*fabs(tmp),0.666666666667);
        if (iwr == 1)if(debug==1)fprintf(myfile," x0 %e\n", x0*HBARC);
        if (IN_PART[2] == 1)
        {
            x0 = sqrt(x0+QUAD(PARTICLE[2].m));
        }
        else
        {
            if (IN_PART[3] == 1)
            {
                x0 = sqrt(x0+QUAD(PARTICLE[3].m));
            }
            else
            {
                x0 = sqrt(x0+QUAD(PARTICLE[4].m));
            }
        }
        if (tmp < 0.) x0 *= -1;
        if (iwr == 1)if(debug==1)fprintf(myfile," x0 %e\n", x0*HBARC);

        f0 = get_n_lep(x0,fc,0);
        ic0 += 1;
        if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %e %e\n", x0*HBARC, f0);

        dx = 0.1;
        if (f0 > 0.)
        {
            xm = x0;
            fm = f0;
            do
            {
                xp = xm;
                fp = fm;
                xm -= dx;

                fm = get_n_lep(xm,fc,0);
                ic0 += 1;
                if (iwr == 1)if(debug==1)fprintf(myfile," xm fm %e %e\n", xm*HBARC, fm);
            }
            while (fm > 0.);
        }
        else
        {
            xp = x0;
            fp = f0;
            do
            {
                xm = xp;
                fm = fp;
                xp += dx;

                fp = get_n_lep(xp,fc,0);
                ic0 += 1;
                if (iwr == 1)if(debug==1)fprintf(myfile," xp fp %e %e\n", xp*HBARC, fp);
            }
            while (fp < 0.);
        }

        if (iwr == 1)
        {
            if(debug==1)fprintf(myfile,"\n xm fm %e %e\n", xm*HBARC, fm);
            if(debug==1)fprintf(myfile," xp fp %e %e\n\n", xp*HBARC, fp);
        }

        tmp = 1.e-09*N_B*V_WS;

        x0old = xp;
        do
        {

            /*if (iwr == 1)if(debug==1)fprintf(myfile,"\n xm fm %f %e\n", xm*HBARC, fm);
            if (iwr == 1)if(debug==1)fprintf(myfile," xp fp %f %e\n\n", xp*HBARC, fp);*/

            x0 = 0.5*(xm+xp);
            /*x0 = (fp*xm-fm*xp)/(fp-fm);*/

            f0 = get_n_lep(x0,fc,0);
            ic += 1;
            if (iwr == 1)if(debug==1)fprintf(myfile," ic x0 f0 %i %f %e\n", ic, x0*HBARC, f0);
            a = f0;
            b = 0.5*(fp-fm);
            c = 0.5*(fp-2.*f0+fm);
            /*printf(" a b c %e %e %e\n", a, b, c);*/
            tmp1 = QUAD(b);
            tmp2 = 4.*a*c;
            tmp3 = 0.5*(xp-xm);
            if (tmp1 > tmp2)
            {
                if (fabs(tmp2/tmp1) > 1.e-06)
                {
                    tmp1 -= tmp2;
                    tmp1 = sqrt(tmp1);
                    x1 = (-b+tmp1)/(2.*c);
                    x2 = (-b-tmp1)/(2.*c);
                }
                else
                {
                    x1 = x2 = -a/b;
                }
                /*printf(" x1 x2 %f %f\n", x1, x2);*/
                if (fabs(x2) < 1.) x1 = x2;
                if (fabs(x1) < 1.)
                {
                    if (fabs(x1) < 0.25)
                    {
                        x2 = x0+2.*x1*tmp3;
                    }
                    else
                    {
                        x2 = x0;
                    }
                    if (fabs(x1) > 0.75)
                    {
                        if (x1 > 0.) x2 = x0+(2.*x1-1)*tmp3;
                        if (x1 < 0.) x2 = x0+(2.*x1+1)*tmp3;
                    }
                    x1 = x0+x1*tmp3;

                    if (f0 < 0.)
                    {
                        xm = x0;
                        fm = f0;
                    }
                    else
                    {
                        xp = x0;
                        fp = f0;
                    }
                    if ((xm < x1) && (x1 < xp))
                    {
                        x0 = x1;
                        f0 = get_n_lep(x0,fc,0);
                        ic += 1;
                        if (iwr == 1)if(debug==1)fprintf(myfile," ic x1 f1 %i %f %e\n", ic, x0*HBARC, f0);
                        if (f0 < 0.)
                        {
                            xm = x0;
                            fm = f0;
                        }
                        else
                        {
                            xp = x0;
                            fp = f0;
                        }
                    }
                    if ((xm < x2) && (x2 < xp))
                    {
                        x0 = x2;
                        f0 = get_n_lep(x0,fc,0);
                        ic += 1;
                        if (iwr == 1)if(debug==1)fprintf(myfile," ic x2 f2 %i %f %e\n", ic, x0*HBARC, f0);
                        if (f0 < 0.)
                        {
                            xm = x0;
                            fm = f0;
                        }
                        else
                        {
                            xp = x0;
                            fp = f0;
                        }
                    }
                }
            }


            if (iwr == 1)if(debug==1)fprintf(myfile," x0 f0 %f %e\n", x0*HBARC, f0);

            if ((xm < x0) && (x0 < xp))
            {
                if (f0 < 0.)
                {
                    xm = x0;
                    fm = f0;
                }
                else
                {
                    xp = x0;
                    fp = f0;
                }
            }

            dx = fabs(x0-x0old);
            x0old = x0;

            if (iwr == 1)if(debug==1)fprintf(myfile," dx xm xp fm fp %e %f %f %e %e\n",
                                                     dx, xm*HBARC, xp*HBARC, fm, fp);

        }
        while (((fabs(f0) > tmp) || (dx > 0.5e-09)) && (ic < 200));

        f0 = get_n_lep(x0,fc,1);
        ic += 1;

    }
    else
    {
        x0 = 0.;
    }

    return x0;
}
/*****************************************************************************/
double get_n_lep(double mu,double fc,int ic)
{
    double fm;
    int il,ip;

    fm = -fc;
    for (il=0; il<3; il++)
    {
        ip = 2+il;
        if (IN_PART[ip] == 1)
        {
            fm += get_rho_part(mu,ip,ic);
        }
    }

    return fm;
}
/*****************************************************************************/
double get_simint(int order_max,double prec,int ifunc,
                  double x_min,double x_max,double rpar,int ipar)
{
    double dx,x,res,tmp1,tmp2,tmp3,err,tmpx;
    int io_max,io,ip_max,ip;

    /* 2010/04/09: more precision needed for convergence */
    prec *= 0.001;

    dx = x_max-x_min;
    res = f_func(x_min,ifunc,rpar,ipar)
          +f_func(x_max,ifunc,rpar,ipar);

    if (order_max > 0)
    {
        io_max = order_max+1;
        io_max = 11;
        io = 0;
        ip_max = 1;
        ip = 1;
        x = x_min+0.5*dx;
        tmp1 = f_func(x,ifunc,rpar,ipar);
        tmp2 = res+4.*tmp1;
        tmp3 = tmp2*dx/6.;
        err = fabs(tmp3);
        for (io=1; io<io_max; io++)
        {
            if ((io < 4) || (err > prec*tmp3))
            {
                tmp3 = tmp2;
                res = res+2.*tmp1;
                dx = 0.5*dx;
                ip_max = 2*ip_max;
                tmp1 = 0.;
                for (ip=0; ip<ip_max; ip++)
                {
                    x = x_min+(0.5+(double)ip)*dx;
                    tmpx = f_func(x,ifunc,rpar,ipar);
                    tmp1 += tmpx;
                }
                tmp2 = res+4.*tmp1;
                err = fabs(dx*(tmp3-0.5*tmp2)/3.);
            }
        }
        res = tmp2*dx/6.;
    }
    else
    {
        res = 0.5*dx*res;
    }

    return res;
}
/*****************************************************************************/
double f_func(double k,int ifunc,double mu,int ipa)
{
    double tmp,k2,tmp1,tmp2;

    switch(ifunc)
    {
    case 1:
    {
        tmp = QUAD(k)*func_fermi(k,mu,ipa);
        break;
    }
    case 2:
    {
        /*tmp = QUAD(k)*func_fermi_p(k,mu,ipa);*/
        break;
    }
    case 3:
    {
        k2 = QUAD(k);
        if (F_MD2 > 0.)
        {
            tmp = k2*fabs(F_MD)/sqrt(k2+F_MD2);
            tmp *= func_fermi(k,mu,ipa);
        }
        else
        {
            tmp = 0.;
        }
        break;
    }
    case 4:
    {
        k2 = QUAD(k);
        tmp = k2*sqrt(k2+F_MD2);
        tmp *= func_fermi(k,mu,ipa);
        break;
    }
    case 5:
    {
        k2 = QUAD(k);
        if (F_MD2 > 0.)
        {
            tmp = k2*k2/sqrt(k2+F_MD2);
            tmp *= func_fermi(k,mu,ipa);
        }
        else
        {
            tmp = k2*k;
        }
        break;
    }
    case 6:
    {
        tmp1 = func_fermi(k,mu,ipa);
        if (tmp1 > 0.)
        {
            tmp2 = tmp1*log(tmp1);
            tmp1 = 1.-tmp1;
            if (tmp1 > 0.)
            {
                tmp2 += tmp1*log(tmp1);
            }
        }
        else
        {
            tmp2 = 0.;
        }
        tmp = QUAD(k)*tmp2;
        break;
    }
    case 11:
    {
        /*
        if (k > 0.) {
          tmp = QUAD(k)*func_bose(k,mu,ipa);
        }
        else {
          tmp = 0.;
          if ((IBX == 1) && (ipa == 0)) tmp = 2.*fabs(F_MD)*TT;
          if ((IBX == -1) && (ipa == 1)) tmp = 2.*fabs(F_MD)*TT;
        }
        */
        k2 = QUAD(k);
        tmp1 = sqrt(k2+F_MD2);
        tmp2 = F_V-mu;
        if (ipa == 0)
        {
            if ((tmp1+tmp2)/TT > 1.e-07)
            {
                tmp = k2*func_bose(k,mu,ipa);
            }
            else
            {
                tmp = 2.*fabs(F_MD)*TT;
            }
        }
        else
        {
            if ((tmp1-tmp2)/TT > 1.e-07)
            {
                tmp = k2*func_bose(k,mu,ipa);
            }
            else
            {
                tmp = 2.*fabs(F_MD)*TT;
            }
        }
        /*if(debug==1)fprintf(FILE_PLOT," %e %e\n", k, tmp);*/
        break;
    }
    case 12:
    {
        /*tmp = QUAD(k)*func_bose_p(k,mu,ipa);*/
        break;
    }
    case 13:
    {
        /*
        tmp = 0.;
        switch (IBX) {
        case 1: {
          if (ipa == 0) tmp = 2.*F_MD*TT;
          break;
        }
        case -1: {
          if (ipa == 1) tmp = 2.*F_MD*TT;
          break;
        }
        default: {
          k2 = QUAD(k);
          tmp = k2*F_MD/sqrt(k2+F_MD2);
          tmp *= func_bose(k,mu,ipa);
        }
        }
        */
        k2 = QUAD(k);
        tmp1 = sqrt(k2+F_MD2);
        tmp2 = F_V-mu;
        if (ipa == 0)
        {
            if ((tmp1+tmp2)/TT > 1.e-07)
            {
                tmp = k2*fabs(F_MD)*func_bose(k,mu,ipa)/tmp1;
            }
            else
            {
                tmp = 2.*fabs(F_MD)*TT;
            }
        }
        else
        {
            if ((tmp1-tmp2)/TT > 1.e-07)
            {
                tmp = k2*fabs(F_MD)*func_bose(k,mu,ipa)/tmp1;
            }
            else
            {
                tmp = 2.*fabs(F_MD)*TT;
            }
        }
        break;
    }
    case 14:
    {
        /*
        if ((IBX != 0) && (k == 0.)) {
          tmp = 2.*F_MD2*TT;
        }
        else {
          k2 = QUAD(k);
          tmp = k2*sqrt(k2+F_MD2);
          tmp *= func_bose(k,mu,ipa);
        }
        */
        k2 = QUAD(k);
        tmp1 = sqrt(k2+F_MD2);
        tmp2 = F_V-mu;
        if (ipa == 0)
        {
            if ((tmp1+tmp2)/TT > 1.e-07)
            {
                tmp = k2*tmp1*func_bose(k,mu,ipa);
            }
            else
            {
                tmp = 2.*F_MD2*TT;
            }
        }
        else
        {
            if ((tmp1-tmp2)/TT > 1.e-07)
            {
                tmp = k2*tmp1*func_bose(k,mu,ipa);
            }
            else
            {
                tmp = 2.*F_MD2*TT;
            }
        }
        /*if(debug==1)fprintf(FILE_PLOT," %e %e\n", k, tmp);*/
        break;
    }
    case 15:
    {
        /*
        if ((IBX != 0) && (k == 0.)) {
          tmp = 0.;
        }
        else {
          k2 = QUAD(k);
          if (F_MD == 0.) {
        tmp = k;
          }
          else {
        tmp = k2/sqrt(k2+F_MD2);
          }
          tmp *= k2*func_bose(k,mu,ipa);
        }
        */
        /* case MD = 0 ? */
        k2 = QUAD(k);
        tmp1 = sqrt(k2+F_MD2);
        tmp2 = F_V-mu;
        if (ipa == 0)
        {
            if ((tmp1+tmp2)/TT > 1.e-07)
            {
                tmp = k2*func_bose(k,mu,ipa)/tmp1;
            }
            else
            {
                tmp = 2.*TT;
            }
        }
        else
        {
            if ((tmp1-tmp2)/TT > 1.e-07)
            {
                tmp = k2*func_bose(k,mu,ipa)/tmp1;
            }
            else
            {
                tmp = 2.*TT;
            }
        }
        break;
    }
    case 16:
    {
        /*
        if ((IBX != 0) && (k == 0.)) {
          tmp = 0.;
        }
        else {
          tmp1 = func_bose(k,mu,ipa);
          if (tmp1 > 0.) {
        tmp2 = tmp1*log(tmp1);
        tmp1 += 1.;
        tmp = QUAD(k)*(tmp2-tmp1*log(tmp1));
          }
          else {
        tmp = 0.;
          }
        }
        */
        k2 = QUAD(k);
        tmp1 = sqrt(k2+F_MD2);
        tmp2 = F_V-mu;
        tmp = 0.;
        if (ipa == 0)
        {
            if ((tmp1+tmp2)/TT > 1.e-07)
            {
                tmp1 = func_bose(k,mu,ipa);
                if (tmp1 > 0.)
                {
                    tmp2 = tmp1*log(tmp1);
                    tmp1 += 1.;
                    tmp = QUAD(k)*(tmp2-tmp1*log(tmp1));
                }
            }
        }
        else
        {
            if ((tmp1-tmp2)/TT > 1.e-07)
            {
                tmp1 = func_bose(k,mu,ipa);
                if (tmp1 > 0.)
                {
                    tmp2 = tmp1*log(tmp1);
                    tmp1 += 1.;
                    tmp = QUAD(k)*(tmp2-tmp1*log(tmp1));
                }
            }
        }
        break;
    }
    default :
    {
        tmp = 1.;
    }
    }

    return tmp;
}
/***************************************************************************/
double get_k0_red(double x)
{
    double tmp,q,i0,t2;

    /* K_0(x)*exp(x)*sqrt(2*x/pi) */

    if (x > 2.)
    {
        q = 2./x;
        tmp = 1.25331414
              +q*(-0.07832358
                  +q*(0.02189568
                      +q*(-0.01062446
                          +q*(0.00587872
                              +q*(-0.00251540
                                  +q*0.00053208)))));
        tmp /= 1.25331414;
    }
    else
    {
        q = 0.25*QUAD(x);
        tmp = -0.57721566
              +q*(0.42278420
                  +q*(0.23069756
                      +q*(0.03488590
                          +q*(0.00262698
                              +q*(0.00010750
                                  +q*0.00000740)))));
        t2 = QUAD(x/3.75);
        i0 = 1.+t2*(3.5156229
                    +t2*(3.0899424
                         +t2*(1.2067492
                              +t2*(0.2659732
                                   +t2*(0.0360768
                                        +t2*0.0045813)))));
        tmp -= i0*log(0.5*x);
        tmp *= sqrt(0.5*PI/x);
        tmp *= exp(x);
    }
    return tmp;
}
/***************************************************************************/
double get_k1_red(double x)
{
    double tmp,q,i1,t2;

    /* K_1(x)*exp(x)*sqrt(2*x/pi) */

    if (x > 2.)
    {
        q = 2./x;
        tmp = 1.25331414
              +q*(0.23498619
                  +q*(-0.03655620
                      +q*(0.01504268
                          +q*(-0.00780353
                              +q*(0.00325614
                                  -q*0.00068245)))));
        tmp /= 1.25331414;
    }
    else
    {
        q = 0.25*QUAD(x);
        tmp = (1.
               +q*(0.15443144
                   +q*(-0.67278579
                       +q*(-0.18156897
                           +q*(-0.01919402
                               +q*(-0.00110404
                                   -q*0.00004686))))))/x;
        t2 = QUAD(x/3.75);
        i1 = (0.5+t2*(0.87890594
                      +t2*(0.51498869
                           +t2*(0.15084934
                                +t2*(0.02658733
                                     +t2*(0.00301532
                                          +t2*0.00032411))))))*x;
        tmp += i1*log(0.5*x);
        tmp *= sqrt(0.5*PI/x);
        tmp *= exp(x);
    }
    return tmp;
}
/***************************************************************************/
double get_k2_red(double x)
{
    double tmp;

    /* K_2(x)*exp(x)*sqrt(2*x/pi) */

    if (x > 0.)
    {
        tmp = get_k0_red(x)+2.*get_k1_red(x)/x;
    }
    else
    {
        tmp = 0.;
    }

    return tmp;
}
/*****************************************************************************/
double func_fermi(double k,double mu,int ipa)
{
    double tmp;

    if (ipa == 1)
    {
        tmp = (F_V-sqrt(QUAD(k)+F_MD2)-mu)/TT;
        if (tmp < 30.)
        {
            tmp = exp(tmp);
            tmp = tmp/(tmp+1.);
        }
        else
        {
            tmp = 1.;
        }
    }
    else
    {
        tmp = (F_V+sqrt(QUAD(k)+F_MD2)-mu)/TT;
        if (tmp < 30.)
        {
            tmp = 1./(exp(tmp)+1.);
        }
        else
        {
            tmp = exp(-tmp);
        }
    }

    return tmp;
}
/*****************************************************************************/
double func_bose(double k,double mu,int ipa)
{
    double tmp;

    if (ipa == 1)
    {
        tmp = (F_V-sqrt(QUAD(k)+F_MD2)-mu)/TT;
        if (tmp > 0.)
        {
            if(debug==1)fprintf(myfile," %e %e %e %e %e %e\n",
                                    F_V*HBARC, F_MD2*QUAD(HBARC), mu*HBARC, k, TT*HBARC, tmp);
            call_error(13);
        }
        tmp = exp(tmp);
        tmp = tmp/(1.-tmp);
    }
    else
    {
        tmp = (F_V+sqrt(QUAD(k)+F_MD2)-mu)/TT;
        if (tmp < 0.)
        {
            if(debug==1)fprintf(myfile," %e %e %e %e %e %e\n",
                                    F_V*HBARC, F_MD2*QUAD(HBARC), mu*HBARC, k, TT*HBARC, tmp);
            call_error(14);
        }
        if (tmp < 30.)
        {
            tmp = 1./(exp(tmp)-1.);
        }
        else
        {
            tmp = exp(-tmp);
        }
    }

    /*
    if ((IBX != 0) && (k == 0.)) {
      tmp = 1.;
    }
    else {
      if (ipa == 1) {
        tmp = (F_V-sqrt(QUAD(k)+F_MD2)-mu)/TT;
        if (tmp > 0.) {
         if(debug==1)fprintf(myfile," %e %e %e %e %e %e\n",
           F_V*HBARC, F_MD2*QUAD(HBARC), mu*HBARC, k, TT*HBARC, tmp);
          call_error(13);
        }
        tmp = exp(tmp);
        tmp = tmp/(1.-tmp);
      }
      else {
        tmp = (F_V+sqrt(QUAD(k)+F_MD2)-mu)/TT;
        if (tmp < 0.) {
         if(debug==1)fprintf(myfile," %e %e %e %e %e %e\n",
           F_V*HBARC, F_MD2*QUAD(HBARC), mu*HBARC, k, TT*HBARC, tmp);
          call_error(14);
        }
        if (tmp < 100.) {
          tmp = 1./(exp(tmp)-1.);
        }
        else {
          tmp = exp(-tmp);
        }
      }
    }
    */

    return tmp;
}
/*****************************************************************************/
int get_simint_bounds_f(double mu,int ipa)
{
    double em,ep,tmp;

    em = disp_rel(0.,ipa);
    tmp = 24.*TT;
    if (ipa == 1)
    {
        ep = mu-tmp;
        if (ep > (em-tmp)) ep = em-tmp;
        if (mu <= (-em))
        {
            em = mu+tmp;
            F_KM = inv_disp_rel(em,ipa);
        }
        else
        {
            F_KM = 0.;
        }
    }
    else
    {
        ep = mu+tmp;
        if (ep < (em+tmp)) ep = em+tmp;
        if (mu >= em)
        {
            em = mu-tmp;
            F_KM = inv_disp_rel(em,ipa);
        }
        else
        {
            F_KM = 0.;
        }
    }
    F_KP = inv_disp_rel(ep,ipa);

    return 0;
}
/*****************************************************************************/
int get_simint_bounds_b(double mu,int ipa)
{
    double em,ep,tmp;

    F_KM = 0.;
    em = disp_rel(0.,ipa);
    tmp = 24.*TT;
    if (ipa == 1)
    {
        ep = em-tmp;
    }
    else
    {
        ep = em+tmp;
    }
    F_KP = inv_disp_rel(ep,ipa);

    return 0;
}
/*****************************************************************************/
double inv_disp_rel(double e,int ipa)
{
    double ml,k,tmp;

    ml = e-F_V;
    k = 0.;
    if (ipa == 1)
    {
        if (ml < 0.)   /* antiparticle */
        {
            tmp = QUAD(ml)-F_MD2;
            if (tmp > 0.) k = sqrt(tmp);
        }
    }
    else   /* particle */
    {
        if (ml > 0.)
        {
            tmp = QUAD(ml)-F_MD2;
            if (tmp > 0.) k = sqrt(tmp);
        }
    }

    return k;
}
/*****************************************************************************/
double disp_rel(double k,int ipa)
{
    double tmp;

    if (ipa == 1)   /* anti-particle */
    {
        tmp = F_V-sqrt(QUAD(k)+F_MD2);
    }
    else   /* particle */
    {
        tmp = F_V+sqrt(QUAD(k)+F_MD2);
    }

    return tmp;
}
/*****************************************************************************/
int get_sources(void)
{
    int ir,imf,ip,ifold,iq,i1;

    /* calculation of source densities */

    /* isoscalar and isovector nucleon densities */
    for (ir=0; ir<NR; ir++)
    {
        /* vector densities */
        DENS_N[0][ir] = DENS[0][ir][2];
        DENS_N[1][ir] = DENS[1][ir][2];
        /* isoscalar */
        DENS_N[2][ir] = DENS_N[0][ir]+DENS_N[1][ir];
        /* isovector */
        DENS_N[3][ir] = DENS_N[0][ir]-DENS_N[1][ir];
        /* scalar densities */
        DENSS_N[0][ir] = DENSS[0][ir][2];
        DENSS_N[1][ir] = DENSS[1][ir][2];
        /* isoscalar */
        DENSS_N[2][ir] = DENSS_N[0][ir]+DENSS_N[1][ir];
        /* isovector */
        DENSS_N[3][ir] = DENSS_N[0][ir]-DENSS_N[1][ir];
    }

    /* RHO_MF and RHOC_MF(condensate) trennen ! */

    /* source densities for photon, omega, sigma, rho fields */
    /* from free nucleons */
    if (0 == 1)
    {
        /* folding with nucleon density */
        ip = 0;
        for (ir=0; ir<NR; ir++)
            DENS_FOLD[ir] = DENS_N[0][ir];
        if (IPHASE > 0) fold_dens(ip);
        for (ir=0; ir<NR; ir++) RHO_MF[0][ir] = DENS_FOLD[ir];
        for (ir=0; ir<NR; ir++)
            DENS_FOLD[ir] = DENS_N[2][ir];
        if (IPHASE > 0) fold_dens(ip);
        for (ir=0; ir<NR; ir++) RHO_MF[1][ir] = DENS_FOLD[ir];
        for (ir=0; ir<NR; ir++)
            DENS_FOLD[ir] = DENSS_N[2][ir];
        if (IPHASE > 0) fold_dens(ip);
        for (ir=0; ir<NR; ir++) RHO_MF[2][ir] = DENS_FOLD[ir];
        for (ir=0; ir<NR; ir++)
            DENS_FOLD[ir] = DENS_N[3][ir];
        if (IPHASE > 0) fold_dens(ip);
        for (ir=0; ir<NR; ir++) RHO_MF[3][ir] = DENS_FOLD[ir];
        for (ir=0; ir<NR; ir++)
            DENS_FOLD[ir] = DENSS_N[3][ir];
        if (IPHASE > 0) fold_dens(ip);
        for (ir=0; ir<NR; ir++) RHO_MF[4][ir] = DENS_FOLD[ir];

        for (ir=0; ir<NR; ir++) RHO_MF[5][ir] = 0.;
    }
    else
    {
        /* without folding of nucleon density */
        for (ir=0; ir<NR; ir++)
        {
            RHO_MF[0][ir] = DENS_N[0][ir];
            RHO_MF[1][ir] = DENS_N[2][ir];
            RHO_MF[2][ir] = DENSS_N[2][ir];
            RHO_MF[3][ir] = DENS_N[3][ir];
            RHO_MF[4][ir] = DENSS_N[3][ir];
            RHO_MF[5][ir] = 0.;
        }
    }

    /* contributions from condensed densities */
    for (imf=0; imf<N_MF; imf++)
    {
        for (ir=0; ir<NR; ir++) RHOC_MF[imf][ir] = 0.;
    }

    /* from leptons */
    for (ip=2; ip<6; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            for (ir=0; ir<NR; ir++)
            {
                RHO_MF[0][ir] -= (DENS[ip][ir][2]);
            }
        }
    }

    /* from table of nuclei, only for homogeneous calculation */
    /* nicht hier !!! */
    if (0 == 1)
    {
        if (IN_TON == 1)
        {
            for (ip=0; ip<NP_NUC; ip++)
            {
                iq = ISORT[ip];
                /*printf(" %i %i %i %i %e\n",
                ip, iq, NUC[iq].a, NUC[iq].z, DENS_TON[iq]);*/
                RHO_MF[0][0] += DENS_TON[iq]*(double)NUC[iq].z;
                RHO_MF[1][0] += DENS_TON[iq]*(double)NUC[iq].a;
                RHO_MF[2][0] += DENSS_TON[iq]*(double)NUC[iq].a;
                RHO_MF[3][0] += DENS_TON[iq]*(double)(NUC[iq].z-NUC[iq].n);
                RHO_MF[4][0] += DENSS_TON[iq]*(double)(NUC[iq].z-NUC[iq].n);
            }
            /*
            if (NP_NUC > 0) {
             if(debug==1)fprintf(myfile," rearrangement contributions! nicht hier!\n");
              exit(0);
            }
            */
        }
    }

    /* from other particles */
    for (ip=9; ip<N_PART; ip++)
    {
        if (IN_PART[ip] == 1)
        {
            /*
            if ((IPHASE > 0) && (PARTICLE[ip].rms > 0.)) {
            ifold = 1;
            }
            else {
            ifold = 0;
            }
                 */
            ifold = 0;
            /* vector density */
            for (ir=0; ir<NR; ir++)
            {
                DENS_FOLD[ir] = DENS[ip][ir][2];
            }
            if (ifold > 0) fold_dens(ip);
            if (PARTICLE[ip].z != 0)
            {
                for (ir=0; ir<NR; ir++)
                {
                    RHO_MF[0][ir] += DENS_FOLD[ir]*(double)PARTICLE[ip].z;
                }
            }
            if (PARTICLE[ip].fo != 0.)
            {
                for (ir=0; ir<NR; ir++)
                    RHO_MF[1][ir] += DENS_FOLD[ir]*PARTICLE[ip].fo;
            }
            if (PARTICLE[ip].fr != 0.)
            {
                for (ir=0; ir<NR; ir++)
                    RHO_MF[3][ir] += DENS_FOLD[ir]*PARTICLE[ip].fr;
            }
            if (PARTICLE[ip].fp != 0.)
            {
                for (ir=0; ir<NR; ir++)
                    RHO_MF[5][ir] += DENS_FOLD[ir]*PARTICLE[ip].fp;
            }
            /* scalar density */
            for (ir=0; ir<NR; ir++)
            {
                DENS_FOLD[ir] = DENSS[ip][ir][2];
            }
            if (ifold > 0) fold_dens(ip);
            if (PARTICLE[ip].fs != 0.)
            {
                for (ir=0; ir<NR; ir++)
                    RHO_MF[2][ir] += DENS_FOLD[ir]*PARTICLE[ip].fs;
            }
            if (PARTICLE[ip].fd != 0.)
            {
                for (ir=0; ir<NR; ir++)
                    RHO_MF[4][ir] += DENS_FOLD[ir]*PARTICLE[ip].fd;
            }

            /* condensate contribution for bosons */
            if (PARTICLE[ip].st < 1)
            {
                for (ir=0; ir<NR; ir++)
                {
                    DENS_FOLD[ir] = DENSC[ip][ir][2];
                }
                if (ifold > 0) fold_dens(ip);
                if (PARTICLE[ip].z != 0)
                {
                    for (ir=0; ir<NR; ir++)
                    {
                        RHOC_MF[0][ir] += DENS_FOLD[ir]*(double)PARTICLE[ip].z;
                    }
                }
                if (PARTICLE[ip].fo != 0.)
                {
                    for (ir=0; ir<NR; ir++)
                        RHOC_MF[1][ir] += DENS_FOLD[ir]*PARTICLE[ip].fo;
                }
                if (PARTICLE[ip].fr != 0.)
                {
                    for (ir=0; ir<NR; ir++)
                        RHOC_MF[3][ir] += DENS_FOLD[ir]*PARTICLE[ip].fr;
                }
                if (PARTICLE[ip].fp != 0.)
                {
                    for (ir=0; ir<NR; ir++)
                        RHOC_MF[5][ir] += DENS_FOLD[ir]*PARTICLE[ip].fp;
                }
                for (ir=0; ir<NR; ir++)
                {
                    DENS_FOLD[ir] = DENSCS[ip][ir][2];
                }
                if (ifold > 0) fold_dens(ip);
                if (PARTICLE[ip].fs != 0.)
                {
                    for (ir=0; ir<NR; ir++)
                        RHOC_MF[2][ir] += DENS_FOLD[ir]*PARTICLE[ip].fs;
                }
                if (PARTICLE[ip].fd != 0.)
                {
                    for (ir=0; ir<NR; ir++)
                        RHOC_MF[4][ir] += DENS_FOLD[ir]*PARTICLE[ip].fd;
                }
            }

        }
    }

    /* densities for density dependence of couplings */
    /* nucleons */
    for (ir=0; ir<NR; ir++)
    {
        RHO_TOT[ir]  = DENS_N[2][ir];
        RHOI_TOT[ir] = DENS_N[3][ir];
    }

    if (DBDEP != 0)
    {
        /* clusters, without condensate contribution so far */
        if (IN_CL == 1)
        {
            for (i1=0; i1<N_CL; i1++)
            {
                ip = 9+i1;
                if (IN_PART[ip] == 1)
                {
                    for (ir=0; ir<NR; ir++)
                    {
                        RHO_TOT[ir]  += DENS[ip][ir][2]*(double)PARTICLE[ip].a;
                        RHOI_TOT[ir] += DENS[ip][ir][2]
                                        *(double)(PARTICLE[ip].z-PARTICLE[ip].n);
                    }
                }
            }
        }
        /* table of nuclei */
        if (IN_TON == 1)
        {
            for (ip=0; ip<NP_NUC; ip++)
            {
                iq = ISORT[ip];
                RHO_TOT[ir]  += DENS_TON[iq]*(double)NUC[iq].a;
                RHOI_TOT[ir] += DENS_TON[iq]*(double)(NUC[iq].z-NUC[iq].n);
            }
        }
    }

    return 0;
}
/*****************************************************************************/
int init_acc_mf(int imf)
{
    int ip,ir,iq;

    if (REC_ACC != (DIM_LA+1)) call_error(300);

    M_ACC_MF = 0;
    for (ip=0; ip<DIM_LA; ip++)
    {
        for (iq=0; iq<DIM_LA; iq++)
        {
            AMF[ip][iq] = 0.;
        }
    }

    iq = 0;
    for (ir=0; ir<NR; ir++)
    {
        VMF[iq][DIM_LA] = MF[imf][ir];
        /*printf(" ip iq %i %i\n", ip, iq);*/
        iq += 1;
    }
    N_ACC_MF = iq;
    /*printf(" N_ACC = %i\n", N_ACC);*/

    for (iq=0; iq<N_ACC_MF; iq++)
    {
        for (ip=0; ip<REC_ACC; ip++)
        {
            FMF[iq][ip] = DFMF[iq][ip] = DVMF[iq][ip] =
                                             UMF[iq][ip] = GUMF[iq][ip] = 0.;
        }
    }

    return 0;
}
/*****************************************************************************/
int get_acc_mf(int imf)
{
    int ip,iq,ir,idx[REC_ACC],ik,in,itmp;
    double relax,tmpf,tmpm;

    for (ip=0; ip<REC_ACC; ip++)
    {
        INDEX_MF[ip] = idx[ip] = (ip+M_ACC)%REC_ACC;
    }

    /*
    for (ip=0; ip<REC_ACC; ip++) {
      if ((idx[ip] < 0) || (idx[ip] > (REC_ACC+1))) {
       if(debug==1)fprintf(myfile," ip idx %i %i\n", idx[ip]);
        call_error(3333);
      }
    }
    */

    M_ACC_MF += 1;

    iq = 0;
    for (ir=0; ir<NR; ir++)
    {
        VMF[iq][idx[0]] = MF[imf][ir];
        iq += 1;
    }

    for (iq=0; iq<N_ACC_MF; iq++)
    {
        FMF[iq][idx[DIM_LA]] = VMF[iq][idx[0]]-VMF[iq][idx[DIM_LA]];
    }

    /* modified Broyden's method */
    /* A. Baran et al., Phys. Rev. C 78 (2008) 014318 */
    relax = 0.6;
    if (M_ACC_MF > 1)
    {
        itmp = DIM_LA-1;
        tmpf = tmpm = 0.;
        for (iq=0; iq<N_ACC_MF; iq++)
        {
            DFMF[iq][idx[itmp]] = FMF[iq][idx[itmp+1]]-FMF[iq][idx[itmp]];
            DVMF[iq][idx[itmp]] = VMF[iq][idx[itmp+1]]-VMF[iq][idx[itmp]];
            tmpm += fabs(VMF[iq][idx[0]]);
            tmpf += QUAD(DFMF[iq][idx[itmp]]);
        }
        /*printf(" tmpm tmpf %e %e\n", tmpm, tmpf);*/
        if (tmpf > 0.) tmpf = 1./sqrt(tmpf);
        for (iq=0; iq<N_ACC_MF; iq++)
        {
            DFMF[iq][idx[itmp]] *= tmpf;
            DVMF[iq][idx[itmp]] *= tmpf;
            UMF[iq][idx[itmp]] = relax*DFMF[iq][idx[itmp]]+DVMF[iq][idx[itmp]];
        }

        for (ik=0; ik<DIM_LA; ik++)
        {
            for (in=0; in<DIM_LA; in++)
            {
                MAT[ik][in] = 0.;
                for (iq=0; iq<N_ACC_MF; iq++)
                    MAT[ik][in] += DFMF[iq][idx[in]]*DFMF[iq][idx[ik]];
            }
        }
        for (ik=0; ik<DIM_LA; ik++) MAT[ik][ik] += 0.0001;

        get_mat_inv(DIM_LA);

        for (ik=0; ik<DIM_LA; ik++)
        {
            CMF[ik] = 0.;
            for (iq=0; iq<N_ACC_MF; iq++)
                CMF[ik] += DFMF[iq][idx[ik]]*FMF[iq][idx[DIM_LA]];
        }

        for (in=0; in<DIM_LA; in++)
        {
            GMF[in] = 0.;
            for (ik=0; ik<DIM_LA; ik++)
                GMF[in] += CMF[ik]*MAT[ik][in];
        }

        for (iq=0; iq<N_ACC_MF; iq++)
        {
            for (ip=0; ip<DIM_LA; ip++)
                GUMF[iq][idx[ip]] = UMF[iq][idx[ip]]*GMF[ip];
        }

    }

    for (iq=0; iq<N_ACC_MF; iq++)
    {
        VMF[iq][idx[0]] = VMF[iq][idx[DIM_LA]]+relax*FMF[iq][idx[DIM_LA]];
    }

    for (iq=0; iq<N_ACC_MF; iq++)
    {
        for (ip=0; ip<DIM_LA; ip++)
            VMF[iq][idx[0]] -= GUMF[iq][idx[ip]];
    }

    iq = 0;
    for (ir=0; ir<NR; ir++)
    {
        MF[imf][ir] = VMF[iq][idx[0]];
        iq += 1;
    }

    return 0;
}
/*****************************************************************************/
int get_mf(int ic)
{
    double tmpo,tmpr,a,c,tmp1,tmp2,n_n,n_p,
           s0[N_MF],s[N_MF][DIM_R],tmpv[DIM_R],tmpn_n[DIM_R],tmpn_p[DIM_R],
           tmpf1[DIM_R],tmpf2[DIM_R],mf2[N_MF][DIM_R],dmf[DIM_R];
    int ir,imf,ip,ngs,i1,i2,iq,iq2,iq3;

    /* calculation of source densities */

    get_sources();

    tmpo = 0.5*LA_OMEGA;
    tmpr = 0.5*LA_RHO;

    i2 = 0;
    do
    {
        i2 += 1;

        if (IN_CL == 1)
        {
            for (imf=0; imf<N_MF; imf++)
            {
                for (ir=0; ir<NR; ir++)
                {
                    mf2[imf][ir] = MF[imf][ir];
                }
            }
        }

        /* pseudo densities */
        /*if (ic > -10) {*/
        if (DBDEP == 0)
        {
            if (ic == 0)
            {
                for (ir=0; ir<NR; ir++)
                {
                    RHO_PS[1][ir] = RHO_TOT[ir];
                    RHO_PS[3][ir] = RHOI_TOT[ir];
                }
            }
            else
            {
                for (ir=0; ir<NR; ir++)
                {
                    RHO_PS[1][ir] = LA_OMEGA*MF[1][ir];
                    RHO_PS[3][ir] = LA_RHO*MF[3][ir];
                }
            }
        }
        else
        {
            for (ir=0; ir<NR; ir++)
            {
                RHO_PS[1][ir] = RHO_TOT[ir];
                RHO_PS[3][ir] = RHOI_TOT[ir];
            }
        }

        /* calculation of new meson fields */

        if (IPHASE > 0)
        {

            /* source functions */
            /* regular contributions */
            for (ir=0; ir<NR; ir++)
            {
                /* photon */
                MF_OLD[0][ir] = MF[0][ir];
                s[0][ir] = G_G*(RHO_MF[0][ir]+RHOC_MF[0][ir]);
                /* mesons */
                /*tmp1 = DENS_N[2][ir];*/
                tmp1 = RHO_TOT[ir];
                get_cpl(tmp1,0);
                for (imf=1; imf<N_MF; imf++)
                {
                    s[imf][ir] = CPL[imf][0]*(RHO_MF[imf][ir]+RHOC_MF[imf][ir]);
                    MF_OLD[imf][ir] = MF[imf][ir];
                }
            }

            /* cluster shift contributions */
            if (1 == 1)
            {
                if ((IN_CL == 1) && (DBDEP == 0))
                {
                    for (ir=0; ir<NR; ir++)
                    {
                        tmpn_n[ir] = 0.5*(RHO_PS[1][ir]-RHO_PS[3][ir]);
                        tmpn_p[ir] = 0.5*(RHO_PS[1][ir]+RHO_PS[3][ir]);
                    }
                    for (i1=0; i1<N_CL; i1++)
                    {
                        ip = 9+i1;
                        if (IN_PART[ip] == 1)
                        {
                            /*
                            for (ir=0; ir<NR; ir++) {
                              DENS_FOLD[ir] = DENSS[ip][ir][2]+DENSCS[ip][ir][2];
                            }
                            fold_dens(ip);
                            for (ir=0; ir<NR; ir++) {
                              tmpf1[ir] = DENS_FOLD[ir];
                            }
                            */
                            for (ir=0; ir<NR; ir++)
                            {
                                tmpf1[ir] = DENSS[ip][ir][2]+DENSCS[ip][ir][2];
                                tmpf2[ir] = DENS[ip][ir][2]+DENSC[ip][ir][2];
                            }
                            for (ir=0; ir<NR; ir++)
                            {
                                n_n = tmpn_n[ir];
                                n_p = tmpn_p[ir];
                                get_be_cl3(i1,n_p,n_n);
                                tmp1 = tmpf1[ir]-XXX*tmpf2[ir];
                                s[1][ir] += tmpo*(DBE_CL[0][i1]+DBE_CL[1][i1])*tmp1;
                                s[3][ir] += tmpr*(DBE_CL[0][i1]-DBE_CL[1][i1])*tmp1;
                            }
                        }
                    }
                }
            }

            /* total source densities */
            for (imf=0; imf<N_MF; imf++)
            {
                for (ir=0; ir<NR; ir++) S_MF[imf][ir] = s[imf][ir];
            }

            /* subtraction of homogeneous background */
            for (imf=0; imf<N_MF; imf++)
            {
                s0[imf] = s[imf][NRP];
                for (ir=0; ir<NR; ir++) s[imf][ir] -= s0[imf];
            }

            for (imf=0; imf<N_MF; imf++)
            {
                ngs = 10;
                if (ic > ngs) ngs = ic;
                /* Greens function */
                for (ir=0; ir<NR; ir++) tmpv[ir] = s[imf][ir];
                get_mf_greens(imf,tmpv);
                /* Gauss-Seidel */
                for (ir=0; ir<NR; ir++) s[imf][ir] *= RVEC[1][ir];
                for (i1=0; i1<ngs; i1++)
                {
                    for (ir=1; ir<NRP; ir++)
                    {
                        MF[imf][ir] = (s[imf][ir+1]*GS_SP[ir]
                                       +s[imf][ir-1]*GS_SM[ir]
                                       +s[imf][ir]*GS_S0[ir]
                                       -MF[imf][ir+1]*GS_AP[imf][ir]
                                       -MF[imf][ir-1]*GS_AM[imf][ir])
                                      /GS_A0[imf][ir];
                    }
                    for (ir=(NR-2); ir>0; ir--)
                    {
                        MF[imf][ir] = (s[imf][ir+1]*GS_SP[ir]
                                       +s[imf][ir-1]*GS_SM[ir]
                                       +s[imf][ir]*GS_S0[ir]
                                       -MF[imf][ir+1]*GS_AP[imf][ir]
                                       -MF[imf][ir-1]*GS_AM[imf][ir])
                                      /GS_A0[imf][ir];
                    }
                }
                /* iteration for improved convergence */
                if (1 == 0)
                {
                    iq = 0;
                    do
                    {
                        iq += 1;
                        for (ir=0; ir<NR; ir++) dmf[ir] = MF[imf][ir];
                        for (ir=1; ir<NRP; ir++)
                        {
                            MF[imf][ir] = (s[imf][ir+1]*GS_SP[ir]
                                           +s[imf][ir-1]*GS_SM[ir]
                                           +s[imf][ir]*GS_S0[ir]
                                           -MF[imf][ir+1]*GS_AP[imf][ir]
                                           -MF[imf][ir-1]*GS_AM[imf][ir])
                                          /GS_A0[imf][ir];
                        }
                        for (ir=(NR-2); ir>0; ir--)
                        {
                            MF[imf][ir] = (s[imf][ir+1]*GS_SP[ir]
                                           +s[imf][ir-1]*GS_SM[ir]
                                           +s[imf][ir]*GS_S0[ir]
                                           -MF[imf][ir+1]*GS_AP[imf][ir]
                                           -MF[imf][ir-1]*GS_AM[imf][ir])
                                          /GS_A0[imf][ir];
                        }
                        tmp1 = 0.;
                        for (ir=0; ir<NR; ir++)
                        {
                            dmf[ir] -= MF[imf][ir];
                            tmp1 += fabs(dmf[ir]);
                            /*printf(" %e %e %e\n", RVEC[1][ir], MF[imf][ir], dmf[ir]);*/
                        }
                        /*printf(" iq tmp1 %i %e\n", iq, tmp1);*/
                    }
                    while ((tmp1 > 1.e-10) && (iq < 200));
                    /*if (iq > 1000) {
                     if(debug==1)fprintf(myfile," imf iq %i %i\n", imf, iq);
                      for (ir=0; ir<NR; ir++) {
                       if(debug==1)fprintf(myfile," %f %e\n", RVEC[1][ir], MF[imf][ir]);
                      }
                      exit(0);
                      }*/
                    /*for (ir=0; ir<NR; ir++) {
                     if(debug==1)fprintf(myfile," %e %e %e\n", RVEC[1][ir], MF[imf][ir], dmf[ir]);
                      }
                     if(debug==1)fprintf(myfile," tmp1 %e\n", tmp1);*/
                    /*exit(0);*/
                }
                /* accelerated convergence ??? */
                if (1 == 0)
                {
                    init_acc_mf(imf);
                    iq = 0;
                    iq2 = 0;
                    iq3 = 0;
                    do
                    {
                        iq += 1;
                        for (ir=0; ir<NR; ir++) dmf[ir] = MF[imf][ir];
                        for (ir=1; ir<NRP; ir++)
                        {
                            MF[imf][ir] = (s[imf][ir+1]*GS_SP[ir]
                                           +s[imf][ir-1]*GS_SM[ir]
                                           +s[imf][ir]*GS_S0[ir]
                                           -MF[imf][ir+1]*GS_AP[imf][ir]
                                           -MF[imf][ir-1]*GS_AM[imf][ir])
                                          /GS_A0[imf][ir];
                        }
                        for (ir=(NR-2); ir>0; ir--)
                        {
                            MF[imf][ir] = (s[imf][ir+1]*GS_SP[ir]
                                           +s[imf][ir-1]*GS_SM[ir]
                                           +s[imf][ir]*GS_S0[ir]
                                           -MF[imf][ir+1]*GS_AP[imf][ir]
                                           -MF[imf][ir-1]*GS_AM[imf][ir])
                                          /GS_A0[imf][ir];
                        }
                        /*if ((iq2 == 1) && (iq3 == 0)) init_acc_mf(imf);
                          if ((iq2 == 1) && (iq3 == 1)) get_acc_mf(imf);*/
                        get_acc_mf(imf);
                        tmp1 = 0.;
                        for (ir=0; ir<NR; ir++)
                        {
                            dmf[ir] -= MF[imf][ir];
                            tmp1 += fabs(dmf[ir]);
                        }
                        if (iq2 == 1) iq3 = 1;
                        /*if ((tmp1 < 1.e-04) || (iq > 100)) iq2 = 1;*/
                        if (tmp1 < 1.e-06) iq2 = 1;
                        /*printf(" iq tmp1 %i %e\n", iq, tmp1);*/
                    }
                    while (tmp1 > 1.e-10);
                    /*printf(" imf iq %i %i\n", imf, iq);*/
                    /*if (imf == 1) {
                      for (ir=0; ir<NR; ir++) {
                     if(debug==1)fprintf(myfile," %f %e\n", RVEC[1][ir], MF[imf][ir]);
                      }
                      exit(0);
                      }*/
                }
            }

            for (imf=0; imf<N_MF; imf++)
            {
                /* boundary conditions and background contribution */
                if (imf == 0)
                {
                    /* Coulomb */
                    if ((IN_PART[2] == 1) || (IN_PART[3] == 1) || (IN_PART[4] == 1))
                    {
                        /* asymptotics for zero net charge included in R_WS */
                        /* xxxxxxxx */
                        if (1 == 0)
                        {
                            /* old */
                            tmp1 = s0[0]/6.;
                            for  (ir=1; ir<NR; ir++) MF[imf][ir] = MF[imf][ir]/RVEC[1][ir]
                                                                       +tmp1*(RVEC[2][NRP]-RVEC[2][ir]);
                        }
                        else
                        {
                            /* new */
                            for (ir=1; ir<NR; ir++)
                            {
                                MF[0][ir] /= RVEC[1][ir];
                            }

                            a = G_G*RHO_MF[0][NRP];
                            c = -a/6.;
                            a *= (F_COUL1+QUAD(R1)/6.);
                            /*a *= (F_COUL1+F_COUL2+QUAD(R1)/6.);*/

                            for (ir=1; ir<NR; ir++)
                            {
                                tmp1 = RVEC[1][ir];
                                MF[0][ir] += a+c*QUAD(tmp1);
                            }

                            /*
                              if (HALLO == 1) {
                             if(debug==1)fprintf(myfile," a c R1 R2 %e %e %e %e\n", a, c, R1, R2);
                             if(debug==1)fprintf(myfile," s0 %e\n", s0[0]);
                              for (ir=1; ir<NR; ir++) {
                              tmp1 = RVEC[1][ir]-R1;
                              MF[0][ir] *= RVEC[1][ir];
                              if(debug==1)fprintf(FILE_PLOT," %f %e %e %e %e\n",
                              RVEC[1][ir], RHO_MF[0][ir],
                              -s[0][ir], (-S_MF[0][ir]*RVEC[1][ir]),
                              MF[0][ir]);
                              }
                              exit(0);
                              }
                            */
                        }
                    }
                    else
                    {
                        /* asymptotics without charge compensation */
                        tmp1 = s0[0]/6.;
                        tmp2 = GS_FP_M[NR-2]*MF[imf][NR-3]
                               +GS_FP_0[NR-2]*MF[imf][NR-2]
                               +GS_FP_P[NR-2]*MF[imf][NRP];
                        for (ir=1; ir<NR; ir++)
                            MF[imf][ir] = MF[imf][ir]/RVEC[1][ir]-tmp2
                                          +tmp1*(RVEC[2][NRP]-RVEC[2][ir]);
                    }
                }
                else
                {
                    for (ir=1; ir<NR; ir++) MF[imf][ir] = MF[imf][ir]/RVEC[1][ir]
                                                              +s0[imf]/MESON[imf].m2;
                }
                MF[imf][0] = (16.*MF[imf][1]-MF[imf][2])/15.;
            }

        }
        else
        {

            n_n = 0.5*(RHO_PS[1][0]-RHO_PS[3][0]);
            n_p = 0.5*(RHO_PS[1][0]+RHO_PS[3][0]);

            if (IN_TON == 1) get_ton_mass(n_p,n_n);

            /* source functions */
            /* regular contributions */
            /* mesons */
            /*tmp1 = DENS_N[2][0];*/
            tmp1 = RHO_TOT[0];
            get_cpl(tmp1,0);
            s[0][0] = 0.;

            for (imf=1; imf<N_MF; imf++)
            {
                s[imf][0] = CPL[imf][0]*(RHO_MF[imf][0]+RHOC_MF[imf][0]);
                MF_OLD[imf][0] = MF[imf][0];
            }

            if (NL == 1)
            {
                s[1][0] = NL_GO*(RHO_MF[1][0]+RHOC_MF[1][0]);
                s[2][0] = NL_GS*(RHO_MF[2][0]+RHOC_MF[2][0]);
                s[3][0] = NL_GR*(RHO_MF[3][0]+RHOC_MF[3][0]);
                s[4][0] = 0.;
                s[5][0] = 0.;
            }

            /* cluster shift contributions */
            if ((IN_CL == 1) && (DBDEP == 0))
            {
                for (i1=0; i1<N_CL; i1++)
                {
                    ip = 9+i1;
                    if (IN_PART[ip] == 1)
                    {

                        get_be_cl3(i1,n_p,n_n);

                        tmp1 = DENSS[ip][0][2]+DENSCS[ip][0][2];
                        tmp1 -= XXX*(DENS[ip][0][2]+DENSC[ip][0][2]);
                        s[1][0] += tmpo*(DBE_CL[0][i1]+DBE_CL[1][i1])*tmp1;
                        s[3][0] += tmpr*(DBE_CL[0][i1]-DBE_CL[1][i1])*tmp1;


                    }
                }
            }

            /* table of nuclei */
            if (IN_TON == 1)
            {
                for (ip=0; ip<NP_NUC; ip++)
                {
                    iq = ISORT[ip];
                    /* regular contributions */
                    s[1][0] += CPL[1][0]*DENS_TON[iq]*(double)NUC[iq].a;
                    s[2][0] += CPL[2][0]*DENSS_TON[iq]*(double)NUC[iq].a;
                    s[3][0] += CPL[3][0]*DENS_TON[iq]*(double)(NUC[iq].z-NUC[iq].n);
                    s[4][0] += CPL[4][0]*DENSS_TON[iq]*(double)(NUC[iq].z-NUC[iq].n);
                    /* mass shift contributions */
                    s[1][0] += NUC[iq].dmdo*(DENSS_TON[iq]-DENS_TON[iq]);
                    s[3][0] += NUC[iq].dmdr*(DENSS_TON[iq]-DENS_TON[iq]);
                }
            }
            /* total source densities */

            for (imf=0; imf<N_MF; imf++) S_MF[imf][0] = s[imf][0];
            MF[0][0] = 0.;
            if (NL == 1)
            {
                MF[1][0] = s[1][0]/(MESON[1].m2+NL_C3*QUAD(MF[1][0]));
                MF[2][0] = s[2][0]/(MESON[2].m2-NL_G2*MF[2][0]+NL_G3*QUAD(MF[2][0]));
                MF[3][0] = s[3][0]/MESON[3].m2;
                MF[4][0] = 0.;
                MF[5][0] = 0.;
            }
            else
            {
                for (imf=1; imf<N_MF; imf++)
                {
                    MF[imf][0] = s[imf][0]/MESON[imf].m2;
                }
            }

        }

        tmp1 = tmp2 = 0.;
        if (IN_CL == 1)
        {

            if (i2 > 2)
            {
                for (imf=0; imf<N_MF; imf++)
                {
                    for (ir=0; ir<NR; ir++)
                    {
                        MF[imf][ir] = 0.5*(MF[imf][ir]+mf2[imf][ir]);
                    }
                }
            }

            for (imf=0; imf<N_MF; imf++)
            {
                for (ir=0; ir<NR; ir++)
                {
                    if (MF[imf][ir] != 0.)
                    {
                        tmp1 += fabs(MF_OLD[imf][ir]/MF[imf][ir]-1.);
                        tmp2 += 1.;
                    }
                }
            }
            if (tmp2 > 0.) tmp1 /= tmp2;
            /*printf(" tmp1 = %e\n", tmp1);*/
            if (i2 > 200)
            {
                if(debug==1)fprintf(myfile," get_mf tmp1 = %e\n", tmp1);
                call_error(40);
            }
        }

    }
    while (tmp1 > 1.e-12);

    /*
    if (NR > 1) {
    for (ir=0; ir<NR; ir++) {
      if(debug==1)fprintf(FILE_PLOT2," %e %e %e %e %e\n",
        RVEC[1][ir], MF[0][ir], MF[1][ir], MF[2][ir], MF[3][ir]);
    }
    exit(0);
    }
    */


    /*printf("i2 tmp1 %i %e\n", i2, tmp1);*/

    /*printf(" MF %e %e %e %e\n", MF[0][0], MF[1][0], MF[2][0], MF[3][0]);*/

    /* test of mean fields*/
    if ((1 == 0) && (IPHASE > 0))
    {
        for (ir=1; ir<(NR-1); ir++)
        {
            if(debug==1)fprintf(FILE_PLOT2," %e %e %e\n",
                                    RVEC[1][ir],
                                    (2.*MF[3][ir]*RVEC[1][ir]
                                     -(MF[3][ir-1]*RVEC[1][ir-1]+MF[3][ir+1]*RVEC[1][ir+1]))
                                    /RVEC[2][1]+MESON[3].m2*MF[3][ir]*RVEC[1][ir],
                                    RVEC[1][ir]*S_MF[3][ir]);
        }
        //  exit(0);

    }

    return 0;
}
/*****************************************************************************/
int get_self_energies(int ic)
{
    double tmp1,tmp_mf[N_MF][DIM_R],tmp_r[DIM_R],n_n,n_p;
    int ir,imf,i1,ip;

    /*if (IN_TON == 1) get_ton_mass(n_p,n_n);*/

    /* calculation of meson fields */

    get_mf(ic);

    /*printf(" MF %e %e %e\n", MF[1][0], MF[2][0], MF[3][0]);*/

    /* calculation of new self-energies */

    if (IPHASE > 0)
    {
        /* inhomogeneous calculation */
        for (ir=0; ir<NR; ir++)
        {
            /*tmp1 = DENS_N[2][ir];*/
            tmp1 = RHO_TOT[ir];
            get_cpl(tmp1,0);
            tmp_mf[0][ir] = G_G*MF[0][ir];
            for (imf=1; imf<N_MF; imf++) tmp_mf[imf][ir] = CPL[imf][0]*MF[imf][ir];
            /* rearrangement contribution from couplings */
            tmp_r[ir] = CPL[1][1]*MF[1][ir]*(RHO_MF[1][ir]+RHOC_MF[1][ir])
                        -CPL[2][1]*MF[2][ir]*(RHO_MF[2][ir]+RHOC_MF[2][ir])
                        +CPL[3][1]*MF[3][ir]*(RHO_MF[3][ir]+RHOC_MF[3][ir])
                        -CPL[4][1]*MF[4][ir]*(RHO_MF[4][ir]+RHOC_MF[4][ir]);
            /*tmp1 = CPL[1][1]*MF[1][ir]*(RHO_MF[1][ir]+RHOC_MF[1][ir])
            -CPL[2][1]*MF[2][ir]*(RHO_MF[2][ir]+RHOC_MF[2][ir])
            +CPL[3][1]*MF[3][ir]*(RHO_MF[3][ir]+RHOC_MF[3][ir])
            -CPL[4][1]*MF[4][ir]*(RHO_MF[4][ir]+RHOC_MF[4][ir]);*/
            /* baryon self-energies */
            VV[0][ir] = VV[1][ir] = tmp_mf[1][ir]+tmp_r[ir];
            /*VV[0][ir] = VV[1][ir] = tmp_mf[1][ir]+tmp1;*/
            VV[0][ir] += tmp_mf[3][ir]+tmp_mf[0][ir];
            VV[1][ir] -= tmp_mf[3][ir];
            SS[0][ir] = SS[1][ir] = tmp_mf[2][ir];
            SS[0][ir] += tmp_mf[4][ir];
            SS[1][ir] -= tmp_mf[4][ir];
            /* hyperon self-energies */
            /*
            if (IN_P[9] == 1) {
            for (i1=0; i1<N_HYP; i1++) {
              VV_H[i1][ir] = tmp_mf[0][ir]*(double)HYPERON[i1].z
                +tmp_mf[1][ir]*HYPERON[i1].fo
                +tmp_mf[3][ir]*HYPERON[i1].fr;
              SS_H[i1][ir] = tmp_mf[2][ir]*HYPERON[i1].fs
                +tmp_mf[4][ir]*HYPERON[i1].fd;
            }
                 }
                 */
            /* thermal meson self-energies */
        }
        /* electron self-energies */
        if (IN_PART[2] == 1)
        {
            for (ir=0; ir<NR; ir++) VV[2][ir] = -G_G*MF[0][ir];
        }
        /* muon self-energies */
        if (IN_PART[3] == 1)
        {
            for (ir=0; ir<NR; ir++) VV[3][ir] = -G_G*MF[0][ir];
        }
        /* tauon self-energies */
        if (IN_PART[4] == 1)
        {
            for (ir=0; ir<NR; ir++) VV[4][ir] = -G_G*MF[0][ir];
        }
        /* cluster self-energies and correction to baryon self-energies */
        if (IN_CL == 1)
        {
            for (i1=0; i1<N_CL; i1++)
            {
                ip = 9+i1;
                if (IN_PART[ip] == 1)
                {
                    for (ir=0; ir<NR; ir++)
                    {
                        n_n = 0.5*(RHO_PS[1][ir]-RHO_PS[3][ir]);
                        n_p = 0.5*(RHO_PS[1][ir]+RHO_PS[3][ir]);
                        get_be_cl3(i1,n_p,n_n);
                        VV[ip][ir] = -BE_CL[i1]*XXX+G_G*MF[0][ir]*(double)PARTICLE[ip].z;
                        SS[ip][ir] = -BE_CL[i1];
                    }
                    /* rearrangement contribution from couplings */
                    if (DBDEP != 0)
                    {
                        for (ir=0; ir<NR; ir++)
                        {
                            VV[ip][ir] += tmp_r[ir]*(double)PARTICLE[ip].a;
                        }
                    }
                    if (PARTICLE[ip].fo != 0.)
                    {
                        for (ir=0; ir<NR; ir++)
                            VV[ip][ir] += tmp_mf[1][ir]*PARTICLE[ip].fo;
                    }
                    if (PARTICLE[ip].fr != 0.)
                    {
                        for (ir=0; ir<NR; ir++)
                            VV[ip][ir] += tmp_mf[3][ir]*PARTICLE[ip].fr;
                    }
                    if (PARTICLE[ip].fs != 0.)
                    {
                        for (ir=0; ir<NR; ir++)
                            SS[ip][ir] += tmp_mf[2][ir]*PARTICLE[ip].fs;;
                    }
                    if (PARTICLE[ip].fd != 0.)
                    {
                        for (ir=0; ir<NR; ir++)
                            SS[ip][ir] += tmp_mf[4][ir]*PARTICLE[ip].fd;
                    }
                }
            }
        }
    }
    else
    {
        /* homogeneous calculation */
        if (NL == 1)
        {
            tmp_mf[1][0] = NL_GO*MF[1][0];
            tmp_mf[2][0] = NL_GS*MF[2][0];
            tmp_mf[3][0] = NL_GR*MF[3][0];
            tmp_mf[4][0] = 0.;
            tmp_mf[5][0] = 0.;
        }
        else
        {
            /*tmp1 = DENS_N[2][0];
            tmp1 = DENS_TOT[0];
            get_cpl(tmp1,0);*/
            for (imf=1; imf<N_MF; imf++)
            {
                tmp_mf[imf][0] = CPL[imf][0]*MF[imf][0];
                /*printf(" %i %e %e\n", imf, CPL[imf][0], MF[imf][0]);*/
            }
        }

        /*printf(" MF %e %e %e\n", MF[1][0], MF[2][0], MF[3][0]);*/

        /* rearrangement contribution from couplings */
        if (NL == 1)
        {
            tmp_r[0] = 0.;
            /*tmp1 = 0.;*/
        }
        else
        {
            tmp_r[0] = CPL[1][1]*MF[1][0]*(RHO_MF[1][0]+RHOC_MF[1][0])
                       -CPL[2][1]*MF[2][0]*(RHO_MF[2][0]+RHOC_MF[2][0])
                       +CPL[3][1]*MF[3][0]*(RHO_MF[3][0]+RHOC_MF[3][0])
                       -CPL[4][1]*MF[4][0]*(RHO_MF[4][0]+RHOC_MF[4][0]);
            /*tmp1 = CPL[1][1]*MF[1][0]*(RHO_MF[1][0]+RHOC_MF[1][0])
            -CPL[2][1]*MF[2][0]*(RHO_MF[2][0]+RHOC_MF[2][0])
            +CPL[3][1]*MF[3][0]*(RHO_MF[3][0]+RHOC_MF[3][0])
            -CPL[4][1]*MF[4][0]*(RHO_MF[4][0]+RHOC_MF[4][0]);*/
        }
        /* baryon self-energies */
        VV[0][0] = VV[1][0] = tmp_mf[1][0]+tmp_r[0];
        /*VV[0][0] = VV[1][0] = tmp_mf[1][0]+tmp1;*/
        VV[0][0] += tmp_mf[3][0];
        VV[1][0] -= tmp_mf[3][0];
        SS[0][0] = SS[1][0] = tmp_mf[2][0];
        SS[0][0] += tmp_mf[4][0];
        SS[1][0] -= tmp_mf[4][0];
        /* 2010/11/24 limitation by zero effective mass */
        if (SS[0][0] > PARTICLE[0].m) SS[0][0] = PARTICLE[0].m-0.5e-05;
        if (SS[1][0] > PARTICLE[1].m) SS[1][0] = PARTICLE[1].m-0.5e-05;
        /* hyperon self-energies */
        if (IN_HYP == 1)
        {
            for (i1=0; i1<N_HYP; i1++)
            {
                ip = 9+N_CL+i1;
                if (IN_PART[ip] == 1)
                {
                    VV[ip][0] = tmp_mf[1][0]*PARTICLE[ip].fo
                                +tmp_mf[3][0]*PARTICLE[ip].fr
                                +tmp_mf[5][0]*PARTICLE[ip].fp;
                    SS[ip][0] = tmp_mf[2][0]*PARTICLE[ip].fs
                                +tmp_mf[4][0]*PARTICLE[ip].fd;
                    /* limitation */
                    if (SS[ip][0] > PARTICLE[ip].m) SS[ip][0] = PARTICLE[ip].m-0.5e-05;
                }
            }
        }
        /* thermal meson self-energies */
        if (IN_TMES == 2)   /* presently not needed */
        {
            for (i1=0; i1<N_TMES; i1++)
            {
                ip = 9+N_CL+N_HYP+i1;
                if (IN_PART[ip] == 1)
                {
                    VV[ip][0] = tmp_mf[1][0]*PARTICLE[ip].fo
                                +tmp_mf[3][0]*PARTICLE[ip].fr
                                +tmp_mf[5][0]*PARTICLE[ip].fp;
                    SS[ip][0] = tmp_mf[2][0]*PARTICLE[ip].fs
                                +tmp_mf[4][0]*PARTICLE[ip].fd;
                    /* limitation */
                    if (SS[ip][0] > PARTICLE[ip].m) SS[ip][0] = PARTICLE[ip].m-0.5e-05;
                }
            }
        }

        /* pseudo densities */
        n_n = 0.5*(RHO_PS[1][0]-RHO_PS[3][0]);
        n_p = 0.5*(RHO_PS[1][0]+RHO_PS[3][0]);

        /* cluster self-energies and correction to baryon self-energies */
        if (IN_CL == 1)
        {
            for (i1=0; i1<N_CL; i1++)
            {
                ip = 9+i1;
                if (IN_PART[ip] == 1)
                {
                    get_be_cl3(i1,n_p,n_n);
                    VV[ip][0] = -BE_CL[i1]*XXX;
                    SS[ip][0] = -BE_CL[i1];
                    /*MM[ip][0] = PARTICLE[ip].m+BE_CL[i1];*/
                    /*printf(" %i %i %f\n", i1, ip, BE_CL[i1]*HBARC);*/
                    /* rearrangement contribution from couplings */
                    if (DBDEP != 0)
                    {
                        VV[ip][0] += tmp_r[0]*(double)PARTICLE[ip].a;
                    }
                    if (PARTICLE[ip].fo != 0.)
                    {
                        VV[ip][0] += tmp_mf[1][0]*PARTICLE[ip].fo;
                    }
                    if (PARTICLE[ip].fr != 0.)
                    {
                        VV[ip][0] += tmp_mf[3][0]*PARTICLE[ip].fr;
                    }
                    if (PARTICLE[ip].fs != 0.)
                    {
                        SS[ip][0] += tmp_mf[2][0]*PARTICLE[ip].fs;
                    }
                    if (PARTICLE[ip].fd != 0.)
                    {
                        SS[ip][0] += tmp_mf[4][0]*PARTICLE[ip].fd;
                    }
                }
            }
        }
        /* self-energies of nuclei */
        if (IN_TON == 1)
        {
            for (i1=0; i1<N_AME11; i1++)
            {
                /* VV_NUC[i1] = 0.; G_G*MF[0][0]*(double)CLUSTER[i1].z;*/
                /* SS_NUC[i1] = 0.; BE_CL[i1];*/
                VV_NUC[i1] = tmp_mf[1][0]*(double)NUC[i1].a;
                VV_NUC[i1] += tmp_mf[3][0]*(double)(NUC[i1].z-NUC[i1].n);
                SS_NUC[i1] = tmp_mf[2][0]*(double)NUC[i1].a;
                SS_NUC[i1] += tmp_mf[4][0]*(double)(NUC[i1].z-NUC[i1].n);
            }
        }
    }

    /*
    if(debug==1)fprintf(myfile," p m S V %e %e %e\n",
     PARTICLE[0].m*HBARC, SS[0][0]*HBARC, VV[0][0]*HBARC);
    if(debug==1)fprintf(myfile," n m S V %e %e %e\n",
     PARTICLE[1].m*HBARC, SS[1][0]*HBARC, VV[1][0]*HBARC);
    if(debug==1)fprintf(myfile," c m S V %e %e %e\n",
     PARTICLE[14].m*HBARC, SS[14][0]*HBARC, VV[14][0]*HBARC);
    */

    /*
    if (NR > 1) {
      for (ir=0; ir<NR; ir++) {
        if(debug==1)fprintf(FILE_PLOT2," %e %e %e %e %e\n",
          RVEC[1][ir], SS[0][ir]*HBARC, SS[1][ir]*HBARC,
          VV[0][ir]*HBARC, VV[1][ir]*HBARC);
      }
      exit(0);
    }
    */

    return 0;
}
/*****************************************************************************/
int init_self_energies(void)
{
    int ir,ip;

    for (ip=0; ip<N_PART; ip++)
    {
        for (ir=0; ir<DIM_R; ir++)
        {
            VV[ip][ir] = SS[ip][ir] = 0.;
            MM[ip][ir] = PARTICLE[ip].m;
        }
    }
    if (IN_TON == 1)
    {
        for (ip=0; ip<N_AME11; ip++) VV_NUC[ip] = SS_NUC[ip] = 0.;
    }

    return 0;
}
/*****************************************************************************/
int init_mf(void)
{
    int ir,imf;

    for (imf=0; imf<N_MF; imf++)
    {
        for (ir=0; ir<DIM_R; ir++) MF[imf][ir] = 0.;
    }

    return 0;
}
/*****************************************************************************/
int init_be_cl(void)
{
    double delta3,tmp1,tmp2,tmp3,eref,sig,c_o,c_s,c_r,c_d,
           fac[2],g_orp,g_orm,g_sdp,g_sdm,facp[2];

    int i1,ip;

    get_cpl(0.,0);

    c_o = QUAD(CPL[1][0])/MESON[1].m2;
    c_s = QUAD(CPL[2][0])/MESON[2].m2;
    c_r = QUAD(CPL[3][0])/MESON[3].m2;
    c_d = QUAD(CPL[4][0])/MESON[4].m2;

    g_orp = c_o+c_r;
    g_sdp = c_s+c_d;
    g_orm = c_o-c_r;
    g_sdm = c_s-c_d;

    facp[0] = -1.5/PARTICLE[0].m;
    facp[1] = -1.5/PARTICLE[1].m;
    fac[0] = 1.+facp[0]*TT;
    fac[1] = 1.+facp[1]*TT;

    TDGGMM[4] = TDGGMM[5] = TDGGMM[6] = TDGGMM[7] = 0.;

    /* pp 1S0 */
    tmp1 = PARTICLE[0].m/sqrt(TPI);
    tmp2 = 0.5*QUAD(PARTICLE[0].g)*CUBE(tmp1);
    GGMM[7] = tmp2*(g_orp-g_sdp*fac[0]*fac[0]);
    if (GGMM[7] != 0.) TDGGMM[7] = -2.*TT*tmp2*g_sdp*facp[0]*fac[0];

    /* nn 1S0 */
    tmp1 = PARTICLE[1].m/sqrt(TPI);
    tmp2 = 0.5*QUAD(PARTICLE[1].g)*CUBE(tmp1);
    GGMM[6] = tmp2*(g_orp-g_sdp*fac[1]*fac[1]);
    if (GGMM[6] != 0.) TDGGMM[6] = -2.*TT*tmp2*g_sdp*facp[1]*fac[1];

    /* pn */
    tmp1 = sqrt(PARTICLE[0].m*PARTICLE[1].m/TPI);
    tmp2 = 0.5*PARTICLE[0].g*PARTICLE[1].g*CUBE(tmp1);
    tmp3 = (facp[0]*fac[1]+fac[0]*facp[1]);
    /* pn 1S0 */
    GGMM[5] = tmp2*(g_orp-g_sdp*fac[0]*fac[1]);
    if (GGMM[5] != 0.) TDGGMM[5] = -TT*tmp2*g_sdp*tmp3;

    /* pn 3S1 */
    GGMM[4] = 2.*tmp2*(g_orm-g_sdm*fac[0]*fac[1])-GGMM[5];
    if (GGMM[4] != 0.) TDGGMM[4] = -TT*2.*tmp2*g_sdm*tmp3-TDGGMM[5];

    for (i1=0; i1<N_CL; i1++)
    {
        delta3 = eref = 0.;
        ip = 9+i1;
        if ((i1 > 3) && (i1 < 8) && (TT > 0.) && (IN_PART[ip] == 1))
        {
            /* 2H and two-body continuum correlations */
            /* continuum contributions with effective range contribution */
            tmp1 = 1./sqrt(2.*CL_MU[i1]*TT*CL_A[i1]);
            tmp2 = 1./sqrt(2.*CL_MU[i1]*TT*CL_B[i1]);
            /* integral */
            VI[i1] = 0.5*CL_C[i1]*experfc(tmp1)/sqrt(CL_A[i1]);
            VI[i1] += 0.5*CL_D[i1]*experfc(tmp2)/sqrt(CL_B[i1]);

            /* T d/dT virial integral*/
            TDVI[i1] = 0.5*CL_C[i1]*tmp1*(1.-sqrt(PI)*tmp1*experfc(tmp1))
                       /(sqrt(PI*CL_A[i1]));
            TDVI[i1] += 0.5*CL_D[i1]*tmp2*(1.-sqrt(PI)*tmp2*experfc(tmp2))
                        /(sqrt(PI*CL_B[i1]));

            if (VI[i1] > 0.)
            {
                sig = 1.;
            }
            else
            {
                sig = -1.;
            }
            if (VI[i1] != 0.)
            {
                eref = -log(sig*VI[i1]);
                delta3 = eref-TDVI[i1]/VI[i1];
                eref *= TT;
                if(debug==1)fprintf(myfile," i1 E_ref %i %e\n", i1, eref*HBARC);
            }
            else
            {
                eref = 1.e+99;
                delta3 = 0.;
            }
            E_RES[i1] = eref;
            DBE_CL[2][i1] = delta3;
        }
        else
        {
            VI[i1] = TDVI[i1] = 0.;
            E_RES[i1] = DBE_CL[2][i1] = 0.;
        }
    }

    return 0;
}
/*****************************************************************************/
int get_be_cl3(int i1,double n_p,double n_n)
{
    double delta,deltap,deltan,deltat,delta4,sig,tmp,tmpt,m,n,n0;
    int ip;

    /* new parametrization of energy shifts in MeV from G. Roepke */

    m = 0.5*(939.565+938.783)/HBARC;
    /*printf("hallo i1 %i\n", i1);*/
    /*TT = 0.01/HBARC;*/

    deltat = 0.;

    ip = 9+i1;
    switch (i1)
    {
    case 1:   /* 3H */
    {
        get_dbe_p_g(i1);
        tmp = DEG[i1];
        tmpt = DEGT[i1];
        /*printf("3H   t %f %f\n", TT*HBARC, DEG[i1]*HBARC);*/
        n0 = PARTICLE[ip].be/tmp;
        n = (2.*n_p+4.*n_n)/3.;
        get_shift(tmp,tmpt,n,n0);
        delta  = SHIFT_F;
        deltap = 2.*SHIFT_FP/3.;
        deltan = 4.*SHIFT_FP/3.;
        deltat = SHIFT_FT;
        break;
    }
    case 2:   /* 3He */
    {
        get_dbe_p_g(i1);
        tmp = DEG[i1];
        tmpt = DEGT[i1];
        /*printf("3He  t %f %f\n", TT*HBARC, DEG[i1]*HBARC);*/
        n0 = PARTICLE[ip].be/tmp;
        n = (4.*n_p+2.*n_n)/3.;
        get_shift(tmp,tmpt,n,n0);
        delta  = SHIFT_F;
        deltap = 4.*SHIFT_FP/3.;
        deltan = 2.*SHIFT_FP/3.;
        deltat = SHIFT_FT;
        break;
    }
    case 3:   /* 4He */
    {
        get_dbe_p_g(i1);
        tmp = DEG[i1];
        tmpt = DEGT[i1];
        n0 = PARTICLE[ip].be/tmp;
        /*printf("4He  T  BE DEG  n0 %f %e %e %e\n",
           TT*HBARC, CLUSTER[i1].be*HBARC, DEG[i1]*HBARC, n0);
           exit(0);*/
        n = n_p+n_n;
        get_shift(tmp,tmpt,n,n0);
        delta  = SHIFT_F;
        deltap = SHIFT_FP;
        deltan = SHIFT_FP;
        deltat = SHIFT_FT;
        break;
    }
    case 8:   /* 3H continuum */
    {
        get_dbe_p_g(i1);
        tmp = DEG[i1];
        tmpt = DEGT[i1];
        /*printf("3H   t %f %f\n", TT*HBARC, DEG[i1]*HBARC);*/
        n0 = PARTICLE[ip].be/tmp;
        n = (2.*n_p+4.*n_n)/3.;
        get_shift(tmp,tmpt,n,n0);
        delta  = SHIFT_F;
        deltap = 2.*SHIFT_FP/3.;
        deltan = 4.*SHIFT_FP/3.;
        deltat = SHIFT_FT;
        break;
    }
    case 9:   /* 3He continuum */
    {
        get_dbe_p_g(i1);
        tmp = DEG[i1];
        tmpt = DEGT[i1];
        /*printf("3He  t %f %f\n", TT*HBARC, DEG[i1]*HBARC);*/
        n0 = PARTICLE[ip].be/tmp;
        n = (4.*n_p+2.*n_n)/3.;
        get_shift(tmp,tmpt,n,n0);
        delta  = SHIFT_F;
        deltap = 4.*SHIFT_FP/3.;
        deltan = 2.*SHIFT_FP/3.;
        deltat = SHIFT_FT;
        break;
    }
    case 10:   /* 4He continuum */
    {
        get_dbe_p_g(i1);
        tmp = DEG[i1];
        tmpt = DEGT[i1];
        n0 = PARTICLE[ip].be/tmp;
        /*printf("4He  T  BE DEG  n0 %f %e %e %e\n",
           TT*HBARC, CLUSTER[i1].be*HBARC, DEG[i1]*HBARC, n0);
           exit(0);*/
        n = n_p+n_n;
        get_shift(tmp,tmpt,n,n0);
        delta  = SHIFT_F;
        deltap = SHIFT_FP;
        deltan = SHIFT_FP;
        deltat = SHIFT_FT;
        break;
    }
    default :   /* 2H and two-body continuum correlations */
    {
        /*get_dbe_p_g(i1);*/
        /*printf("2H   t %f %f\n", TT*HBARC, DEG[i1]*HBARC);*/
        /*TT = 20./HBARC;*/
        get_dbe_p_j();
        tmp = DEJ;
        tmpt = DEJT;
        /*printf("2H   t %f %f\n", TT*HBARC, DEJ*HBARC);*/
        n0 = PARTICLE[9].be/DEJ;
        /*printf(" BE T  n  %f %f %f\n", PARTICLE[9].be*HBARC, TT*HBARC, n0);
          exit(0);*/
        switch (i1)
        {
        case 6:   /* nn 1S0 */
        {
            n = 2.*n_n;
            get_shift(tmp,tmpt,n,n0);
            delta  = SHIFT_F;
            deltap = 0.;
            deltan = 2.*SHIFT_FP;
            deltat = SHIFT_FT;
            break;
        }
        case 7:   /* pp 1S0 */
        {
            n = 2.*n_p;
            get_shift(tmp,tmpt,n,n0);
            delta  = SHIFT_F;
            deltap = 2.*SHIFT_FP;
            deltan = 0.;
            deltat = SHIFT_FT;
            break;
        }
        default :   /* np 3S1 & 1S0 */
        {
            n = n_p+n_n;
            get_shift(tmp,tmpt,n,n0);
            delta  = SHIFT_F;
            deltap = SHIFT_FP;
            deltan = SHIFT_FP;
            deltat = SHIFT_FT;
        }
        }
    }
    }

    /* hello, no medium shifts */
    /*printf(" no medium shifts delta\n");
    if(debug==1)fprintf(myfile,"delta deltap deltan E_RES %f %f %f %f \n",
     delta*HBARC, deltap*HBARC, deltan*HBARC, E_RES[i1]*HBARC);
     delta = deltap = deltan = deltat = 0.;*/

    /* modified degeneracy factor */
    delta4 = 0.;
    if ((i1 > 3) && (i1 < 8) && (TT > 0.) && (IN_PART[ip] == 1))
    {
        if (VI[i1] > 0.)
        {
            sig = 1.;
        }
        else
        {
            sig = -1.;
        }
        if (VI[i1] != 0.)
        {
            if (1 == 1)
            {
                /* with RMF correction */
                tmp = PARTICLE[ip].m+E_RES[i1];
                delta4 = PARTICLE[ip].g0*VI[i1]+sqrt(TT/CUBE(tmp))*GGMM[i1];
            }
            else
            {
                /* without RMF corrections */
                tmp = PARTICLE[ip].m;
                delta4 = PARTICLE[ip].g0*VI[i1];
            }
            /*printf(" i1 GGMM %i %e\n", i1, GGMM[i1]);*/
            PARTICLE[ip].g = delta4/(sig*VI[i1]);
            PARTICLE[ip].tdgdt = TDVI[i1]/VI[i1]-0.5; /*-1.5*DBE_CL[2][i1]/tmp;*/
            if (GGMM[i1] != 0.)
            {
                PARTICLE[ip].tdgdt -= TDGGMM[i1]/GGMM[i1];
            }
            PARTICLE[ip].tdgdt *= sig*PARTICLE[ip].g0/PARTICLE[ip].g-1.;
            /*printf(" g0 g %e %e\n",
            PARTICLE[ip].g0, PARTICLE[ip].g);*/
        }
        else
        {
            PARTICLE[ip].g = 0.;
            /* hallo: Ableitung kann ungleich Null sein */
            PARTICLE[ip].tdgdt = 0.;
        }
    }
    else
    {
        /* no modification */
        /*PARTICLE[ip].g = PARTICLE[ip].g0;*/
        PARTICLE[ip].tdgdt = 0.;
    }

    /*printf(" ip g0 g GGMM E_RES VI %i %e %e %e %e %e\n",
      ip, PARTICLE[ip].g0, PARTICLE[ip].g, GGMM[i1], E_RES[i1]*HBARC, VI[i1]);*/

    /* vacuum binding energy already in mass! */

    /*BE_CL[i1] = -delta-E_RES[i1];
    DBE_CL[0][i1] = -deltap;
    DBE_CL[1][i1] = -deltan;
    DBE_CL[3][i1] = -deltat;*/

    /*delta = deltap = deltan = deltat = 0.;*/
    /*E_RES[i1] = DBE_CL[2][i1] = 0.;*/

    BE_CL[i1] = delta+E_RES[i1];
    DBE_CL[0][i1] = deltap;
    DBE_CL[1][i1] = deltan;
    DBE_CL[3][i1] = deltat;

    /*printf(" m_p m_n m BE %f %f %f %f\n",
     PARTICLE[0].m*HBARC, PARTICLE[1].m*HBARC,
     PARTICLE[ip].m*HBARC, BE_CL[i1]*HBARC);*/

    /*DBE_CL[2][i1] = DBE_CL[3][i1] = 0.;*/
    /*DBE_CL[0][i1] = DBE_CL[1][i1];
      BE_CL[i1] = 0./HBARC;*/

    /*printf(" m_p m_n m BE %f %f %f %f\n",
     PARTICLE[0].m*HBARC, PARTICLE[1].m*HBARC,
     PARTICLE[ip].m*HBARC, BE_CL[i1]*HBARC);*/

    /*printf(" m delta E_res BE %f %f %f %f\n",
     PARTICLE[0].m*HBARC, delta*HBARC,
     E_RES[i1]*HBARC, BE_CL[i1]*HBARC);*/

    /*
    for (ip=0; ip<N_PART; ip++) {
     if(debug==1)fprintf(myfile,"uuu %i %i %i %e\n",
       ip, PARTICLE[ip].a, PARTICLE[ip].z, PARTICLE[ip].mu*HBARC);
    }
    exit(0);
    */

    return 0;
}
/*****************************************************************************/
int get_dbe_p_g(int icl)
{
    double tmp0,tmp1;

    tmp0 = TT+DEPG[1][icl];
    tmp1 = sqrt(tmp0);

    DEG[icl] = DEPG[0][icl]/(tmp0*tmp1);
    DEGT[icl] = -1.5*DEG[icl]/tmp0;

    return 0;
};
/*****************************************************************************/
int get_dbe_p_j(void)
{
    double a,ra,b,b2,tmp;

    if (TT > 0.)
    {
        a = 1.+DEPJ[0]/TT;
        ra = sqrt(a);
        b = DEPJ[1];
        b2 = QUAD(b);
        tmp = 1./ra-sqrt(PI)*b*exp(b2*a)*erfc(b*ra);
        tmp *= DEPJ[2]/(TT*sqrt(TT));
    }
    else
    {
        a = DEPJ[0];
        b = DEPJ[1];
        tmp = 2.*QUAD(b)*a*sqrt(a);
        tmp = DEPJ[2]/tmp;
    }
    DEJ = tmp;

    /* temperature derivative not yet included */
    if (TT > 0.)
    {
        DEJT = -(1.5+DEPJ[0]*b2/TT)*DEJ/TT
               +DEPJ[2]*DEPJ[0]*0.5/(a*ra*pow(TT,3.5));
    }
    else
    {
        /* ??? */
        DEJT = 0.75*DEPJ[2]*(1.5/QUAD(b)-1.)/(QUAD(b)*pow(DEPJ[0],2.5));
    }

    return 0;
}
/****************************************************************************/
int get_shift(double de,double det,double x,double x0)
{
    double arg1,arg2,z,b,c,xs;

    /* particle dependent shift */
    z = 1.;
    /*x_zero = x0*(sqrt(2.*z+1.)-1.)/z;*/
    arg1 = x/x0;
    SHIFT_F  = x0*(1.+0.5*z*arg1)*arg1;
    SHIFT_FP = 1.+z*arg1;
    SHIFT_FT = x*(1.+z*arg1);

    /*arg1 = exp(arg1);
      SHIFT_F  = x0*arg1;
      SHIFT_FP = arg1;
      SHIFT_FT = x*arg1;*/

    if (1 == 0)
    {
        xs = 0.15;
        arg2 = x/xs;
        if (arg2 < 1.)
        {
            SHIFT_F  = x/(1.-arg2);
            SHIFT_FP = 1./QUAD(1.-arg2);
            SHIFT_FT = SHIFT_F;
        }
        else
        {
            SHIFT_F = 1.e+20;
            SHIFT_FP = 0.;
            SHIFT_FT = 0.;
        }
    }

    SHIFT_F *= de;
    SHIFT_FP *= de;
    SHIFT_FT *= det;

    /* density dependent shift */

    /*b = 70.22407666;*/
    /*
      b = 35.11203833;
      c = b*TT;
      arg1 = x/RHOREF;
      arg2 = exp(-0.5*QUAD(arg1));
      SHIFT_F += c*(1.-arg2);
      SHIFT_FP += c*arg2*arg1/RHOREF;
      SHIFT_FT += b*(1.-arg2);
    */
    if (0 == 1)
    {
        b = 3.;
        c = TT*27.63102112; /* ln(10^12)*/
        arg1 = RHOREF/x;
        arg2 = tanh(b*(1.-arg1));
        z = c*(1.+arg2);
        SHIFT_F += z;
        SHIFT_FP += c*(1.-QUAD(arg2))*b*arg1/x;
        SHIFT_FT += z/TT;
    }

    return 0;
}
/*****************************************************************************/
int init_dens(void)
{
    double a,tmp,n_0,z_0,
           tmpv_n[DIM_R],tmpv_p[DIM_R],r_n,r_p,tmp_n,tmp_p;
    int imf,ip,ir,ia;

    /* initialisation of densities */
    for (ir=0; ir<DIM_R; ir++)
    {
        RHO_TOT[ir] = RHOI_TOT[ir] = 0.;
        for (imf=0; imf<N_MF; imf++)
        {
            RHO_MF[imf][ir] = RHOC_MF[imf][ir] = RHO_PS[imf][ir] = 0.;
        }
        for (ip=0; ip<N_PART; ip++)
        {
            for (ia=0; ia<3; ia++)
            {
                DENS[ip][ir][ia] = DENSS[ip][ir][ia] = 0.;
                DENSC[ip][ir][ia] = DENSCS[ip][ir][ia] = 0.;
            }
        }
    }

    /* only free nucleons */
    if (IPHASE > 0)
    {
        /* inhomogeneous calculation */
        a = 1.2; /*0.7;*/
        if (N_B > 0.)
        {
            z_0 = N_B*Y_Q;
            n_0 = N_B-z_0;
        }
        else
        {
            z_0 = n_0 = 0.;
        }
        switch (IPHASE)
        {
        case 2:   /* bubble */
        {
            r_n = 10.;
            tmp = CUBE(r_n)+3.*AA_WS/(FPI*0.15);
            r_p = pow(tmp,0.333333333333);
            /*printf(" AA_WS R_WS r_n r_p %f %f %f %f\n", AA_WS, R_WS, r_n, r_p);*/
            for (ir=0; ir<NR; ir++)
            {
                tmpv_p[ir] = tmpv_n[ir] =
                                 (1.-f_ws(RVEC[1][ir],r_n,a))*f_ws(RVEC[1][ir],r_p,a);
                /*printf(" r f %f %f\n", RVEC[1][ir], tmpv_p[ir]);*/
            }
            /*exit(0);*/
            break;
        }
        case 3:   /* hole */
        {
            tmp = CUBE(R_WS)-3.*AA_WS/(FPI*0.15);
            if (tmp > 0.)
            {
                r_p = r_n = pow(tmp,0.333333333333);
            }
            else
            {
                r_p = r_n = 0.;
            }
            /*printf(" r_n r_p %f %f\n", r_n, r_p);
            exit(0);*/
            for (ir=0; ir<NR; ir++)
            {
                tmpv_p[ir] = 1.-f_ws(RVEC[1][ir],r_p,a);
                tmpv_n[ir] = 1.-f_ws(RVEC[1][ir],r_n,a);
            }
            break;
        }
        default:   /* nucleus */
        {
            if ((z_0 < 0.075) && (ZZ_WS > 0.))
            {
                tmp = 3.*ZZ_WS/(FPI*(0.075-z_0));
                r_p = pow(tmp,0.333333333333);
            }
            else
            {
                r_p = R1;
            }

            if ((n_0 < 0.075) && (NN_WS > 0.))
            {
                tmp = 3.*NN_WS/(FPI*(0.075-n_0));
                r_n = pow(tmp,0.333333333333);
            }
            else
            {
                r_n = R1;
            }
            /* 2011/04/21 */
            /*if (r_n < r_p) {
            r_p = r_n;
            }
            else {
            r_n = r_p;
            }*/
            /*tmp = 0.5*(r_n+r_p);
            r_n = r_p = tmp;*/
            for (ir=0; ir<NR; ir++)
            {
                tmpv_p[ir] = f_ws(RVEC[1][ir],r_p,a);
                tmpv_n[ir] = f_ws(RVEC[1][ir],r_n,a);
            }
        }
        }

        tmp_n = 1./(FPI*r_integrate(1,2,tmpv_n));
        tmp_p = 1./(FPI*r_integrate(1,2,tmpv_p));

        for (ir=0; ir<NR; ir++)
        {
            DENS[0][ir][0] = ZZ_WS*tmp_p*tmpv_p[ir];
            DENS[1][ir][0] = NN_WS*tmp_n*tmpv_n[ir];
        }
        for (ir=0; ir<NR; ir++)
        {
            tmp = DENS[0][ir][0]+DENS[1][ir][0];
            tmp = 1./(1.+4.*QUAD(tmp));
            DENSS[0][ir][0] = tmp*DENS[0][ir][0];
            DENSS[1][ir][0] = tmp*DENS[1][ir][0];
        }

        /* electrons and muons */
        if ((IN_PART[2]  == 1) || (IN_PART[3] == 1))
        {
            tmp = ZZ_WS/V_WS;
            if (IN_PART[2] == 1)
            {
                for (ir=0; ir<NR; ir++) DENS[2][ir][0] = tmp;
            }
            else
            {
                for (ir=0; ir<NR; ir++) DENS[3][ir][0] = tmp;
            }
        }
    }
    else
    {
        /* homogeneous calculation */
        DENS[0][0][0] = ZZ_WS;
        DENS[1][0][0] = NN_WS;
        /* besser skalieren !*/
        tmp = 1./(1.+4.*QUAD(AA_WS));
        for (ip=0; ip<2; ip++)
        {
            DENSS[ip][0][0] = tmp*DENS[ip][0][0];
        }

        /* electrons and muons */
        if ((IN_PART[2]  == 1) || (IN_PART[3] == 1))
        {
            if (IN_PART[2] == 1)
            {
                DENS[2][0][0] = ZZ_WS;
            }
            else
            {
                DENS[3][0][0] = ZZ_WS;
            }
        }

    }

    for (ir=0; ir<NR; ir++)
    {
        for (ip=0; ip<4; ip++)
        {
            DENS[ip][ir][2] = DENS[ip][ir][0];
            DENSS[ip][ir][2] = DENSS[ip][ir][0];
        }
    }

    if (IN_TON == 1)
    {
        NP_NUC = 0;
        for (ip=0; ip<N_NUC; ip++)
        {
            DENS_TON[ip] = DENSS_TON[ip] = 0.;
        }
    }

    return 0;
}
/*****************************************************************************/
int get_mf_greens(int imf,double *tmps)
{
    int ir,irp,ix;
    double m,arg,argp,tmpv[DIM_R];

    /* calculation with Green's functions */
    if (imf == 0)
    {
        /* Coulomb */
        MF[imf][0] = 0.;
        for (ix=1; ix<(NX+1); ix++)
        {
            ir = 2*ix;

            for (irp=0; irp<(ir+1); irp++) tmpv[irp] = RVEC[2][irp]*tmps[irp];
            MF[imf][ir] = simpson_mod(ir,tmpv)/RVEC[1][ir];

            for (irp=ir; irp<NR; irp++) tmpv[irp-ir] = RVEC[1][irp]*tmps[irp];
            MF[imf][ir] += simpson_mod((NR-ir-1),tmpv);
        }
        irp = NRP;
        for (ix=0; ix<(NX+1); ix++)
        {
            ir = 2*ix;
            MF[imf][ir] -= MF[imf][irp];
            MF[imf][ir] *= RVEC[1][ir];
        }
    }
    else
    {
        /* mesons */
        m = MESON[imf].m;
        MF[imf][0] = 0.;
        for (ix=1; ix<(NX+1); ix++)
        {
            ir = 2*ix;
            arg = m*RVEC[1][ir];

            for (irp=0; irp<(ir+1); irp++)
            {
                argp = m*RVEC[1][irp];
                tmpv[irp] = RVEC[1][irp]*tmps[irp]
                            *(exp(-arg+argp)-exp(-arg-argp));
            }
            MF[imf][ir] = simpson_mod(ir,tmpv);

            for (irp=ir; irp<NR; irp++)
            {
                argp = m*RVEC[1][irp];
                tmpv[irp-ir] = RVEC[1][irp]*tmps[irp]
                               *(exp(-argp+arg)-exp(-argp-arg));
            }
            MF[imf][ir] += simpson_mod((NR-ir-1),tmpv);
        }
        for (ix=1; ix<(NX+1); ix++)
        {
            ir = 2*ix;
            MF[imf][ir] /= (2.*m);
        }
    }

    /* interpolation with Gauss-Seidel */
    for (ir=0; ir<NR; ir++) tmps[ir] *= RVEC[1][ir];
    for (ix=0; ix<NX; ix++)
    {
        ir = 2*ix+1;
        MF[imf][ir] = (tmps[ir+1]*GS_SP[ir]
                       +tmps[ir-1]*GS_SM[ir]
                       +tmps[ir]*GS_S0[ir]
                       -MF[imf][ir+1]*GS_AP[imf][ir]
                       -MF[imf][ir-1]*GS_AM[imf][ir])
                      /GS_A0[imf][ir];
    }

    return 0;
}
/*****************************************************************************/
double simpson_mod(int np,double *vec)
{
    int ir, irp, n;
    double tmp1, tmp2=0., tmp3=0.;

    if (np > 1)
    {
        n = np/2;
        for (ir=0; ir<(n-1); ir++)
        {
            irp = 2*ir+1;
            tmp2 += vec[irp];
            tmp3 += vec[irp+1];
        }
        tmp2 += vec[2*n-1];
        tmp1 = 4.*tmp2+2.*tmp3+vec[0]+vec[np];
    }
    else
    {
        tmp1 = 0.;
    }

    return (tmp1*DR/3.);
}
/*****************************************************************************/
int fold_dens(int ip)
{
    double a2,a,r,rp,rr,c,n,nf,tmp,tmpv[DIM_R],tmpvp[DIM_R];
    int ir,irp;

    /* Achtung: aeusserer Beitrag */
    for (ir=0; ir<NR; ir++) tmpv[ir] = DENS_FOLD[ir];
    n = r_integrate(0,2,tmpv);

    if (n != 0.)
    {
        a = PARTICLE[ip].ag;
        if (a > 0.)
        {
            a2 = QUAD(a);
            c = 1./(a2*a*pow(PI,1.5));

            rr = RVEC[1][NRP];
            for (ir=0; ir<NR; ir++)
            {
                /*ir = 400; {*/
                r = RVEC[1][ir];
                for (irp=0; irp<NR; irp++)
                {
                    rp = RVEC[1][irp];
                    tmpvp[irp] = DENS_FOLD[irp]*func_fold(r,rp,a2);
                    /*if(debug==1)fprintf(FILE_PLOT," %e %e\n", rp, tmpvp[irp]);*/
                }
                tmpv[ir] = r_integrate(-1,2,tmpvp);
                if (DENS_FOLD[NRP] != 0.)
                {
                    /* extension assuming constant density */
                    for (irp=0; irp<NR; irp++)
                    {
                        rp = RVEC[1][irp]+rr;
                        tmpvp[irp] = func_fold(r,rp,a2)*QUAD(rp);
                        /*if(debug==1)fprintf(FILE_PLOT," %e %e\n", rp, tmpvp[irp]);*/
                    }
                    tmpv[ir] += DENS_FOLD[NRP]*trapez(tmpvp);
                }
                tmpv[ir] *= c*TPI;
                /*tmpv[ir] = c*TPI*(r_integrate(0,2,tmpvp)
                  +func_fold2(r,rr,a)*DENS_FOLD[NRP]);*/
                /*if(debug==1)fprintf(FILE_PLOT," %e %e\n", r, tmpv[ir]);*/
            }
            nf = r_integrate(0,2,tmpv);
            tmp = n/nf;
            for (ir=0; ir<NR; ir++)
            {
                DENS_FOLD[ir] = tmp*tmpv[ir];
            }
        }
    }

    /*
    if ((n > 0.) && (ip == 12) && (DENS_FOLD[NRP] < 0.)) {
      for (ir=0; ir<NR; ir++) {
        if(debug==1)fprintf(FILE_PLOT," %e %e\n", RVEC[1][ir], DENS_FOLD[ir]);
      }
      exit(0);
    }
    */

    return 0;
}
/*****************************************************************************/
double func_fold(double r,double rp,double a2)
{
    double tmp;

    tmp = r*rp;
    if (tmp > 0.)
    {
        tmp = 0.5*a2*(exp(-QUAD(r-rp)/a2)-exp(-QUAD(r+rp)/a2))/tmp;
    }
    else
    {
        tmp = 2.*exp(-(QUAD(r)+QUAD(rp))/a2);
    }

    return tmp;
}
/*****************************************************************************/
double r_integrate(int ic,int ip,double *tmpv)
{
    double tmp1,tmp2,vec[DIM_R];
    int ir,irp;

    irp = NRP;
    if (ip > 2) call_error(19);
    if (ic < 0)
    {
        for (ir=0; ir<NR; ir++)
        {
            vec[ir] = RVEC[ip][ir]*tmpv[ir];
        }
        tmp1 = trapez(vec);
        tmp2 = 0.;
    }
    else
    {
        if ((ic == 1) && (ip == 2))
        {
            for (ir=0; ir<NR; ir++)
            {
                vec[ir] = RVEC[ip][ir]*(tmpv[ir]-tmpv[irp]);
            }
            tmp1 = trapez(vec);
            tmp2 = tmpv[irp]*pow(R1,(double)(ip+1))/(double)(ip+1);
            if (R2 > R1)
            {
                tmp2 += tmpv[irp]*VX/FPI;
            }
        }
        else
        {
            for (ir=0; ir<NR; ir++)
            {
                vec[ir] = RVEC[ip][ir]*(tmpv[ir]-tmpv[irp]);
            }
            tmp1 = trapez(vec);
            tmp2 = tmpv[irp]*pow(R1,(double)(ip+1))/(double)(ip+1);
        }
    }

    return (tmp1+tmp2);
}
/*****************************************************************************/
double trapez (double *vec)
{
    int ir;
    double tmp;

    tmp = 0.;
    for (ir=0; ir<NR; ir++)
    {
        tmp += vec[ir];
    }
    tmp -= 0.5*(vec[0]+vec[NRP]);

    return (tmp*DR);
}
/*****************************************************************************/
double f_ws(double r,double rr, double a)
{
    double tmp;

    tmp = 1./(1.+exp((r-rr)/a));

    return tmp;
}
/*****************************************************************************/
int fit_radius2(int iph,double r_in)
{
    int iwr;

    iwr = 1;

    /*CI_R = fit_aa4(iph,r_in);*/
    CI_R = fit_aa4_new(iph,r_in);
    if((int)CI_R == 999) return 999;

    if (iwr == 1)if(debug==1)fprintf(myfile," CI_R = %i CI_RMF = %i\n", CI_R, CI_RMF);

    return 0;
}
/*****************************************************************************/
double solve_a(int ic, int iph,double a)
{
    int iwr;
    double r;

    iwr = 1;

    r = pow((3.*a/(FPI*N_B)),0.333333333333);
    /*printf(" a r %e %e\n", a, r);*/
    if (IMFS[iph] == 0) XCL = 0;
    if((int)solve_rmf(iph,r) == 999) return 999;
    IMFS[iph] = 1;
    XCL = 1;
    if((int)solve_rmf(iph,r) == 999) return 999;
    /*IMFS[iph] = 0;*/
    if (iwr == 1)
    {
        if(debug==1)fprintf(myfile," ic   r a f   %i  %f %f %f   %i\n",
                                ic, r, a, F_DENS*HBARC/N_B, ICV_RMF);
        /*if(debug==1)fprintf(FILE_PLOT,
          " %f %f %f\n", r, a, F_DENS*HBARC/N_B);*/
    }

    return F_DENS;
}
/*****************************************************************************/
int fit_aa4_new(int iph,double r)
{
    double ax,fx,aa[20],ff[20],da,ffmin,tmp1,tmp2,icmax;
    int ic,ix,ip,ixmin,idx[20],dim,dimp,i1,ipmax,icv[20];

    icmax = 100;
    ipmax = 2;
    ic = 0;

    if (1 == 0)
    {
        /*if (TT > 0.) {*/
        /* no fit of a */
        /* hallo */
        r = 14.805;

        /* beta equilibrium EoS DD2 original */
        if (N_B < 2.0e-04)
        {
            tmp1 = 0.1*(log(N_B)+20.);
            tmp2 = 7.64024+tmp1*(-3.36463+tmp1*(0.167564+tmp1*(-0.151532+tmp1*0.142839)));
            r = exp(tmp2);
        }
        else
        {
            tmp1 = 0.1*(log(N_B)+10.);
            tmp2 = 74.0789+tmp1*(-188.195+tmp1*(525.464+tmp1*(-823.498+tmp1*420.850)));
            r = tmp2;
        }
        /*printf(" N_B r %e %e\n", N_B, r);
          exit(0);*/

        if(debug==1)fprintf(myfile," *** radius fixed to r = %f fm\n", r);

        /*r = 19.966054;*/

        ax = N_B*FPI*CUBE(r)/3.;
        fx = solve_a(ic,iph,ax);
        ic += 1;
    }
    else
    {
        /* with fit of a */

        /*r = 20.;*/
        /*r = 28.;*/

        dim = 20;
        dimp = dim-1;

        ax = N_B*FPI*CUBE(r)/3.;

        for (i1=0; i1<1; i1++)
        {
            if (i1 == 0)
            {
                da = 10.;
                ax = N_B*FPI*CUBE(r)/3.;
                da = 0.005*ax;
                if (iph == 2) da *= 4.;
                if(debug==1)fprintf(myfile," x da r %f %f %f\n", ax, da, r);
            }
            else
            {
                da = 0.5;
            }

            /* hallo */
            da *= 2.;

            for (ix=0; ix<dim; ix++)
            {
                aa[ix] = ff[ix] = 0.;
                idx[ix] = icv[ix] = 0;
            }

            aa[dimp] = ax;
            ff[dimp] = solve_a(ic,iph,ax);
            if((int)ff[dimp] == 999) return 999;
            idx[dimp] = 1;
            icv[ix] = ICV_RMF;
            ic += 1;

            ip = 0;
            do
            {
                ax = aa[dimp]+da;
                fx = solve_a(ic,iph,ax);
                if((int)fx == 999) return 999;
                ic += 1;
                for (ix=0; ix<dimp; ix++)
                {
                    aa[ix] = aa[ix+1];
                    ff[ix] = ff[ix+1];
                    idx[ix] = idx[ix+1];
                    icv[ix] = icv[ix+1];
                }
                aa[dimp] = ax;
                ff[dimp] = fx;
                idx[dimp] = 1;
                icv[ix] = ICV_RMF;
                if (ip > 0)
                {
                    if ((ff[dimp] > ff[dimp-1]) && (ff[dimp-1] > ff[dimp-2])) ip += 1;
                }
                else
                {
                    if (ff[dimp] > ff[dimp-1]) ip += 1;
                }
                /*printf(" ip = %i\n", ip);*/
            }
            while ((ip < ipmax) && (ic < icmax));
            /*
            for (ix=0; ix<dim; ix++) {
            if(debug==1)fprintf(myfile," %i %f %f %i\n", ix, aa[ix], ff[ix]*HBARC/N_B, idx[ix]);
            }
                 */
            ip = 0;
            ix = dim-1;
            do
            {
                ix -= 1;
                if (idx[ix] == 0)
                {
                    aa[ix] = ax = aa[ix+1]-da;
                    ff[ix] = solve_a(ic,iph,ax);
                    if((int)ff[ix] == 999) return 999;
                    idx[ix] = 1;
                    icv[ix] = ICV_RMF;
                    ic += 1;
                }
                if (ip > 0)
                {
                    if ((ff[ix] > ff[ix+1]) && (ff[ix+1] > ff[ix+2])) ip += 1;
                }
                else
                {
                    if (ff[ix] > ff[ix+1]) ip += 1;
                }
            }
            while ((ip < ipmax) && (ix > 0) && (ic < icmax));
            /*
            for (ix=0; ix<dim; ix++) {
            if(debug==1)fprintf(myfile," %i %f %f %i\n", ix, aa[ix], ff[ix]*HBARC/N_B, idx[ix]);
            }
                 */
            if (ip < ipmax)
            {
                do
                {
                    ax = aa[0]-da;
                    fx = solve_a(ic,iph,ax);
                    if((int)fx == 999) return 999;
                    ic += 1;
                    for (ix=dimp; ix>0; ix--)
                    {
                        aa[ix] = aa[ix-1];
                        ff[ix] = ff[ix-1];
                        idx[ix] = idx[ix-1];
                        icv[ix] = icv[ix-1];
                    }
                    aa[0] = ax;
                    ff[0] = fx;
                    idx[0] = 1;
                    icv[0] = ICV_RMF;
                    if (ip > 0)
                    {
                        if ((ff[0] > ff[1]) && (ff[1] > ff[2])) ip += 1;
                    }
                    else
                    {
                        if (ff[0] > ff[1]) ip += 1;
                    }
                }
                while ((ip < ipmax) && (ic < icmax));
            }

            for (ix=0; ix<dim; ix++)
            {
                if (idx[ix] == 1)
                    if(debug==1)fprintf(myfile," %i %f %f %i\n", ix, aa[ix], ff[ix]*HBARC/N_B, idx[ix]);
            }

            ixmin = -1;
            ffmin = 1.e99;
            for (ix=0; ix<dim; ix++)
            {
                if (ff[ix] > 0.)
                {
                    if (ff[ix] < ffmin)
                    {
                        ffmin = ff[ix];
                        ixmin = ix;
                    }
                }
            }

            if(debug==1)fprintf(myfile," ixmin = %i\n", ixmin);
            if ((ixmin < 2) || (ixmin > (dim-3)))
            {
                ax = aa[ixmin];
            }
            else
            {
                /*
                 if(debug==1)fprintf(myfile," ixmin aa ffmin %i %f %f\n", ixmin, aa[ixmin], ffmin*HBARC/N_B);
                  */

                /*
                fx = (2.*(ff[ixmin+2]-ff[ixmin-2])*(ff[ixmin+1]+ff[ixmin-1])
                  -(ff[ixmin+1]-ff[ixmin-1])*(ff[ixmin+2]+ff[ixmin-2]))
                  /(4.*(ff[ixmin+2]-ff[ixmin-2])-2.*(ff[ixmin+1]-ff[ixmin-1]));

                ax = aa[ixmin]-(ff[ixmin+1]-ff[ixmin-1])*da
                  /(ff[ixmin+1]-2.*fx+ff[ixmin-1]);
                if(debug==1)fprintf(myfile," ax fx %f %f\n", ax, fx*HBARC/N_B);
                */

                tmp1 = (2.*(ff[ixmin+2]-ff[ixmin-2])*(ff[ixmin+1]+ff[ixmin-1])
                        -(ff[ixmin+1]-ff[ixmin-1])*(ff[ixmin+2]+ff[ixmin-2]));
                tmp2 = 4.*(ff[ixmin+2]-ff[ixmin-2])-2.*(ff[ixmin+1]-ff[ixmin-1]);
                tmp1 = (ff[ixmin+1]+ff[ixmin-1])*tmp2-2.*tmp1;

                if (tmp1 != 0.)
                {
                    ax = aa[ixmin]-(ff[ixmin+1]-ff[ixmin-1])*da*tmp2/tmp1;
                    if (ax < aa[ixmin-1]) ax = aa[ixmin];
                    if (ax > aa[ixmin+1]) ax = aa[ixmin];
                }
                else
                {
                    ax = aa[ixmin];
                }
            }
            fx = solve_a(ic,iph,ax);
            if((int)fx == 999) return 999;
            /*printf("\n ax fx %f %f\n", ax, fx*HBARC/N_B);*/
        }
        /*printf(" ic = %i\n", ic);*/
    }
    /*if (ic > (icmax-1)) ICV_RMF = 0;*/

    return ic;
}
/*****************************************************************************/
int fit_aa4(int iph,double r)
{
    double a0,am,ap,f0,fm,fp,da,r0,rp,rm,a,b,c,ax,fx,rx,xp,xm,xp2,xm2,xpm,
           a_err,a0_old;
    int ic,iwr,ia;

    iwr = 1;
    ic = 0;


    r = 14.810954; /* y = 0.5 */
    r = 16.668882; /* y = 0.4, T = 0 MeV */
    r = 15.837687; /* y = 0.4, T = 5 MeV */
    /*if (TT > 0.) r = 1.;*/

    if (1 == 0)
    {
        /* no fit of a */
        r0 = r;
        a0 = N_B*FPI*CUBE(r0)/3.;
    }
    else
    {
        /* with fit of a */

        a_err = 1.e-04;

        r0 = r;
        a0 = N_B*FPI*CUBE(r0)/3.;


        for (ic=0; ic<30; ic++)
        {
            a0 = 180.+1.*(double)ic;
            f0 = solve_a(ic,iph,a0);
        }
        // exit(0);
        return 0;
        f0 = solve_a(ic,iph,a0);
        ic += 1;


        da = 1.;

        ap = a0+da;
        fp = solve_a(ic,iph,ap);
        ic += 1;


        if (fp < f0)
        {
            da = 2.;
            do
            {
                am = a0;
                fm = f0;
                rm = r0;
                a0 = ap;
                f0 = fp;
                r0 = rp;
                ap += da;
                da += 2.*da;
                if (da > 500.) da = 500.;
                fp = solve_a(ic,iph,ap);
                ic += 1;

            }
            while (fp < f0);
        }
        else
        {
            am = a0-da;
            if (am < 0.) am = 0.5*a0;
            fm = solve_a(ic,iph,am);
            ic += 1;

            if (fm < f0)
            {
                do
                {
                    ap = a0;
                    fp = f0;
                    rp = r0;
                    a0 = am;
                    f0 = fm;
                    r0 = rm;
                    if ((am-da) < 1.)
                    {
                        am *= 0.5;
                    }
                    else
                    {
                        am -= da;
                        da += 2.*da;
                        if (da > 250) da = 250;
                    }
                    fm = solve_a(ic,iph,am);
                    ic += 1;

                }
                while ((fm < f0) && (am > 1.));
            }
        }

        /*da = ap-am;*/

        do
        {
            a0_old = a0;
            ia = 0;

            if ((fp > f0) && (fm > f0))
            {

                xp = ap-a0;
                xm = am-a0;
                xp2 = QUAD(xp);
                xm2 = QUAD(xm);
                xpm = xp-xm;

                a = f0;
                b = -((fp*xm2-fm*xp2)/xpm+a*(xp+xm))/(xp*xm);
                c = ((fp*xm-fm*xp)/xpm+a)/(xp*xm);
                /*printf(" a b c %e %e %e\n", a, b, c);*/


                /*b = 0.5*(fp-fm);
                c = 0.5*(fp-2.*f0+fm);
                da = 0.5*(ap-am);*/

                ax = a0-0.5*b/c;

                /*if ((ax > 0.5*(am+a0)) && (ax < 0.5*(a0+ap))) {*/
                if ((ax > am) && (ax < ap))
                {
                    fx = solve_a(ic,iph,ax);
                    ic += 1;

                    if (fx < f0)
                    {
                        if (ax < a0)
                        {
                            ap = a0;
                            fp = f0;
                            rp = r0;
                            a0 = ax;
                            f0 = fx;
                            r0 = rx;
                            ia = 1;
                        }
                        else
                        {
                            if (ax > a0)
                            {
                                am = a0;
                                fm = f0;
                                rm = r0;
                                a0 = ax;
                                f0 = fx;
                                r0 = rx;
                                ia = 1;
                            }
                        }
                    }
                    else
                    {
                        /*if (fx > f0) {*/
                        if (ax < a0)
                        {
                            am = ax;
                            fm = fx;
                            rm = rx;
                        }
                        else
                        {
                            /*if (ax > a0) {*/
                            ap = ax;
                            fp = fx;
                            rp = rx;
                        }
                    }
                }
            }

            ax = 0.5*(ap+a0);
            fx = solve_a(ic,iph,ax);
            ic += 1;

            if (fx < f0)
            {
                am = a0;
                fm = f0;
                rm = r0;
                a0 = ax;
                f0 = fx;
                r0 = rx;
                ia = 1;
            }
            else
            {
                ap = ax;
                fp = fx;
                rp = rx;
            }

            ax = 0.5*(am+a0);
            fx = solve_a(ic,iph,ax);
            ic += 1;

            if (fx < f0)
            {
                ap = a0;
                fp = f0;
                rp = r0;
                a0 = ax;
                f0 = fx;
                r0 = rx;
                ia = 1;
            }
            else
            {
                am = ax;
                fm = fx;
                rm = rx;
            }


            if (ia == 1)
            {
                da = fabs(a0-a0_old);
            }
            else
            {
                da = fabs(ap-am);
            }


            /*} while (da > a_err);*/
            /*} while ((((da/a0) > a_err) || (dmu > 1.e-08)) && (ic < 100));*/
        }
        while (((da/a0) > a_err) && (ic < 100));

        if (ic > 99) ICV_RMF = 0;

        if (ISOL == 0) return ic;

    }

    f0 = solve_a(ic,iph,a0);
    ic += 1;

    return ic;
}
/*****************************************************************************/
int fit_aa4_old(int iph,double r)
{
    double a0,am,ap,f0,fm,fp,da,r0,rp,rm,a,b,c,ax,fx,rx,xp,xm,xp2,xm2,xpm,
           a_err,a0_old,mu_b,mu_q,mu_b_old,mu_q_old,dmu;
    int ic,iwr,ia,iwrp;

    iwr = 1;
    iwrp = 0;
    ic = 0;

    mu_b = mu_q = 0.;

    r = 14.810954; /* y = 0.5 */
    r = 16.668882; /* y = 0.4, T = 0 MeV */
    r = 15.837687; /* y = 0.4, T = 5 MeV */
    /*if (TT > 0.) r = 1.;*/

    if (1 == 0)
    {
        /* no fit of a */
        r0 = r;
        a0 = N_B*FPI*CUBE(r0)/3.;
    }
    else
    {
        /* with fit of a */

        a_err = 1.e-04;

        r0 = r;
        a0 = N_B*FPI*CUBE(r0)/3.;
        solve_rmf(iph,r0);
        ic += 1;
        f0 = F_DENS;
        mu_b = MU_BB;
        mu_q = MU_QQ;
        dmu = 1.;

        if (iwr == 1)if(debug==1)fprintf(myfile," r0 a0 f0   %f %f %f   %f %f %f\n",
                                                 r0, a0, f0*HBARC/N_B,
                                                 mu_b*HBARC, mu_q*HBARC, dmu*HBARC);

        da = 1.;

        ap = a0+da;
        rp = pow((3.*ap/(FPI*N_B)),0.333333333333);
        mu_b_old = mu_b;
        mu_q_old = mu_q;
        solve_rmf(iph,rp);
        ic += 1;
        fp = F_DENS;
        mu_b = MU_BB;
        mu_q = MU_QQ;
        dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
        if (iwr == 1)if(debug==1)fprintf(myfile," rp ap fp   %f %f %f   %f %f %f\n",
                                                 rp, ap, fp*HBARC/N_B,
                                                 mu_b*HBARC, mu_q*HBARC, dmu*HBARC);

        if (fp < f0)
        {
            da = 2.;
            do
            {
                am = a0;
                fm = f0;
                rm = r0;
                a0 = ap;
                f0 = fp;
                r0 = rp;
                ap += da;
                da += 2.*da;
                if (da > 500.) da = 500.;
                rp = pow((3.*ap/(FPI*N_B)),0.333333333333);
                mu_b_old = mu_b;
                mu_q_old = mu_q;
                solve_rmf(iph,rp);
                ic += 1;
                fp = F_DENS;
                mu_b = MU_BB;
                mu_q = MU_QQ;
                dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
                if (iwr == 1)if(debug==1)fprintf(myfile," rp ap fp   %f %f %f   %f %f %f\n",
                                                         rp, ap, fp*HBARC/N_B,
                                                         mu_b*HBARC, mu_q*HBARC, dmu*HBARC);
            }
            while (fp < f0);
        }
        else
        {
            am = a0-da;
            if (am < 0.) am = 0.5*a0;
            rm = pow((3.*am/(FPI*N_B)),0.333333333333);
            mu_b_old = mu_b;
            mu_q_old = mu_q;
            solve_rmf(iph,rm);
            ic += 1;
            fm = F_DENS;
            mu_b = MU_BB;
            mu_q = MU_QQ;
            dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
            if (iwr == 1)if(debug==1)fprintf(myfile," rm am fm   %f %f %f   %f %f %f\n",
                                                     rm, am, fm*HBARC/N_B,
                                                     mu_b*HBARC, mu_q*HBARC, dmu*HBARC);
            if (fm < f0)
            {
                do
                {
                    ap = a0;
                    fp = f0;
                    rp = r0;
                    a0 = am;
                    f0 = fm;
                    r0 = rm;
                    if ((am-da) < 1.)
                    {
                        am *= 0.5;
                    }
                    else
                    {
                        am -= da;
                        da += 2.*da;
                        if (da > 250) da = 250;
                    }
                    rm = pow((3.*am/(FPI*N_B)),0.333333333333);
                    mu_b_old = mu_b;
                    mu_q_old = mu_q;
                    solve_rmf(iph,rm);
                    ic += 1;
                    fm = F_DENS;
                    mu_b = MU_BB;
                    mu_q = MU_QQ;
                    dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
                    if (iwr == 1)if(debug==1)fprintf(myfile," rm ap fm   %f %f %f   %f %f %f\n",
                                                             rm, am, fm*HBARC/N_B,
                                                             mu_b*HBARC, mu_q*HBARC, dmu*HBARC);
                }
                while ((fm < f0) && (am > 1.));
            }
        }

        if (iwrp == 1)
        {
            if(debug==1)fprintf(myfile,"\n rm am fm %f %f %e\n", rm, am, fm*HBARC/N_B);
            if(debug==1)fprintf(myfile," r0 a0 f0 %f %f %e\n", r0, a0, f0*HBARC/N_B);
            if(debug==1)fprintf(myfile," rp ap fp %f %f %e\n\n", rp, ap, fp*HBARC/N_B);
        }


        /*da = ap-am;*/

        do
        {
            a0_old = a0;
            ia = 0;

            if ((fp > f0) && (fm > f0))
            {

                xp = ap-a0;
                xm = am-a0;
                xp2 = QUAD(xp);
                xm2 = QUAD(xm);
                xpm = xp-xm;

                a = f0;
                b = -((fp*xm2-fm*xp2)/xpm+a*(xp+xm))/(xp*xm);
                c = ((fp*xm-fm*xp)/xpm+a)/(xp*xm);
                /*printf(" a b c %e %e %e\n", a, b, c);*/


                /*b = 0.5*(fp-fm);
                c = 0.5*(fp-2.*f0+fm);
                da = 0.5*(ap-am);*/

                ax = a0-0.5*b/c;

                /*if ((ax > 0.5*(am+a0)) && (ax < 0.5*(a0+ap))) {*/
                if ((ax > am) && (ax < ap))
                {
                    rx = pow((3.*ax/(FPI*N_B)),0.333333333333);
                    mu_b_old = mu_b;
                    mu_q_old = mu_q;
                    solve_rmf(iph,rx);
                    ic += 1;
                    fx = F_DENS;
                    mu_b = MU_BB;
                    mu_q = MU_QQ;
                    dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
                    if (iwr == 1)if(debug==1)fprintf(myfile," rx ax fx   %f %f %f   %f %f %f\n",
                                                             rx, ax, fx*HBARC/N_B,
                                                             mu_b*HBARC, mu_q*HBARC, dmu*HBARC);
                    if (fx < f0)
                    {
                        if (ax < a0)
                        {
                            ap = a0;
                            fp = f0;
                            rp = r0;
                            a0 = ax;
                            f0 = fx;
                            r0 = rx;
                            ia = 1;
                        }
                        else
                        {
                            if (ax > a0)
                            {
                                am = a0;
                                fm = f0;
                                rm = r0;
                                a0 = ax;
                                f0 = fx;
                                r0 = rx;
                                ia = 1;
                            }
                        }
                    }
                    else
                    {
                        /*if (fx > f0) {*/
                        if (ax < a0)
                        {
                            am = ax;
                            fm = fx;
                            rm = rx;
                        }
                        else
                        {
                            /*if (ax > a0) {*/
                            ap = ax;
                            fp = fx;
                            rp = rx;
                        }
                    }
                }
            }

            ax = 0.5*(ap+a0);
            rx = pow((3.*ax/(FPI*N_B)),0.333333333333);
            mu_b_old = mu_b;
            mu_q_old = mu_q;
            solve_rmf(iph,rx);
            ic += 1;
            fx = F_DENS;
            mu_b = MU_BB;
            mu_q = MU_QQ;
            dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
            if (iwr == 1)if(debug==1)fprintf(myfile," rx ax fx   %f %f %f   %f %f %f\n",
                                                     rx, ax, fx*HBARC/N_B,
                                                     mu_b*HBARC, mu_q*HBARC, dmu*HBARC);
            if (fx < f0)
            {
                am = a0;
                fm = f0;
                rm = r0;
                a0 = ax;
                f0 = fx;
                r0 = rx;
                ia = 1;
            }
            else
            {
                ap = ax;
                fp = fx;
                rp = rx;
            }

            ax = 0.5*(am+a0);
            rx = pow((3.*ax/(FPI*N_B)),0.333333333333);
            mu_b_old = mu_b;
            mu_q_old = mu_q;
            solve_rmf(iph,rx);
            ic += 1;
            fx = F_DENS;
            mu_b = MU_BB;
            mu_q = MU_QQ;
            dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
            if (iwr == 1)if(debug==1)fprintf(myfile," rx ax fx   %f %f %f   %f %f %f\n",
                                                     rx, ax, fx*HBARC/N_B,
                                                     mu_b*HBARC, mu_q*HBARC, dmu*HBARC);
            if (fx < f0)
            {
                ap = a0;
                fp = f0;
                rp = r0;
                a0 = ax;
                f0 = fx;
                r0 = rx;
                ia = 1;
            }
            else
            {
                am = ax;
                fm = fx;
                rm = rx;
            }


            if (ia == 1)
            {
                da = fabs(a0-a0_old);
            }
            else
            {
                da = fabs(ap-am);
            }

            if (iwrp == 1)
            {
                if(debug==1)fprintf(myfile,"\n rm am fm %f %f %e\n", rm, am, fm*HBARC/N_B);
                if(debug==1)fprintf(myfile," r0 a0 f0 %f %f %e   %f\n", r0, a0, f0*HBARC/N_B, da);
                if(debug==1)fprintf(myfile," rp ap fp %f %f %e\n\n", rp, ap, fp*HBARC/N_B);
            }

            /*} while (da > a_err);*/
        }
        while ((((da/a0) > a_err) || (dmu > 1.e-08)) && (ic < 100));

        if (ic > 99) ICV_RMF = 0;

        if (ISOL == 0) return ic;

    }

    r0 = pow((3.*a0/(FPI*N_B)),0.333333333333);
    mu_b_old = mu_b;
    mu_q_old = mu_q;
    /*printf(" mu_b mu_q %e %e\n", mu_b*HBARC, mu_q*HBARC);*/
    solve_rmf(iph,r0);
    ic += 1;
    f0 = F_DENS;
    mu_b = MU_BB;
    mu_q = MU_QQ;
    dmu = fabs(mu_b-mu_b_old)+fabs(mu_q-mu_q_old);
    if (iwr == 1)if(debug==1)fprintf(myfile," r0 a0 f0   %f %f %f   %f %f %f\n",
                                             r0, a0, f0*HBARC/N_B,
                                             mu_b*HBARC, mu_q*HBARC, dmu*HBARC);

    return ic;
}
/****************************************************************************/
int init_discr_ws_cell(int iph)
{
    double tmp,a,b,c,d,xp,xm,tmp2,r_max;
    int ir,imf,nrp;

    r_max = 20.;
    if (iph > 1) r_max = 40.;
    if (R_WS > r_max)
    {
        NR = DIM_R;
        R1 = r_max;
    }
    else
    {
        R1 = R_WS;
        if (R1 > 20.)
        {
            NR = (int)(R1/0.1);
            NR = (2*(NR/2))+1;
            if (NR > DIM_R) NR = DIM_R;
        }
        else
        {
            NR = 201;
        }
    }
    R2 = R_WS;

    /*printf(" R1 R2 %f %f\n", R1, R2);*/

    V1 = FPI*CUBE(R1)/3.;
    V_WS = V2 = FPI*CUBE(R2)/3.;
    VX = V2-V1;

    /* Coulomb correction */
    if (R2 > R1)
    {
        F_COUL1 = -QUAD(R2-R1)*(1.+2.*R2/R1)/6.;
        tmp = R1/R2;
        tmp2 = QUAD(tmp);
        F_COUL2 = 0.2*(1.-QUAD(tmp2)*tmp)-tmp2+tmp2*tmp;
        F_COUL2 *= -0.5*QUAD(R2)/(1.-tmp*tmp2);
    }
    else
    {
        F_COUL1 = F_COUL2 = 0.;
    }

    if (NR > DIM_R) call_error(3);

    /* radial grid */
    nrp = NR-1;
    DR = R1/(double)nrp;
    for (ir=0; ir<NR; ir++)
    {
        RVEC[0][ir] = 1.;
        RVEC[1][ir] = DR*(double)ir;
        RVEC[2][ir] = QUAD(RVEC[1][ir]);
    }

    /* Gauss-Seidel iteration */
    for (ir=1; ir<nrp; ir++)
    {
        xp = RVEC[1][ir+1]-RVEC[1][ir];
        xm = RVEC[1][ir-1]-RVEC[1][ir];
        a = xp-xm;
        b = 3.*xp*xm-QUAD(xp)-QUAD(xm);
        tmp = QUAD(xp)-QUAD(xm);
        tmp2 = xp*xm;
        c = tmp+tmp2;
        d = tmp-tmp2;
        GS_S0[ir] = a*b;
        GS_SM[ir] = xp*c;
        GS_SP[ir] = xm*d;
        GS_A0[0][ir] = -12.*a;
        GS_AM[0][ir] =  12.*xp;
        GS_AP[0][ir] = -12.*xm;
        for (imf=1; imf<N_MF; imf++)
        {
            GS_A0[imf][ir] = a*(-12.+b*MESON[imf].m2);
            GS_AM[imf][ir] = xp*( 12.+c*MESON[imf].m2);
            GS_AP[imf][ir] = xm*(-12.+d*MESON[imf].m2);
        }
    }

    for (ir=1; ir<nrp; ir++)
    {
        xp = RVEC[1][ir+1]-RVEC[1][ir];
        xm = RVEC[1][ir-1]-RVEC[1][ir];
        a = xp-xm;

        GS_FP_P[ir] = -xm/(xp*a);
        GS_FP_0[ir] = -(xp+xm)/(xp*xm);
        GS_FP_M[ir] =  xp/(xm*a);

        GS_FPP_P[ir] =  2./(xp*a);
        GS_FPP_0[ir] =  2./(xp*xm);
        GS_FPP_M[ir] = -2./(xm*a);
    }

    return 0;
}
/*****************************************************************************/
int get_parameter_tables(int ir_y,int ir_n,int ir_t,
                         double d_y_q,double f_n_b,double t_ref,double f_t,char *cnsfile)
{
    int ni_y_q,i_y_q,ni_n_b,i_n_b,ni_t,i_t;
    double y_q,n_b,t;

    /* Y_q */
    if(debug==1) FILE_PLOT = fopen(changestring("yq",cnsfile,".dat"),"w");
    ni_y_q = 1;
    if(debug==1)fprintf(FILE_PLOT," %i\n", ni_y_q);
    ni_y_q = 100;
    if(debug==1)fprintf(FILE_PLOT," %i\n", ni_y_q);
    for (i_y_q=1; i_y_q<(ni_y_q+1); i_y_q++)
    {
        y_q = d_y_q*(double)i_y_q;
        if(debug==1)fprintf(FILE_PLOT," %14.7e\n", y_q);
    }
    if(debug==1) fclose(FILE_PLOT);
    /* n_b */
    if(debug==1)FILE_PLOT = fopen(changestring("nb",cnsfile,".dat"),"w");
    ni_n_b = 1;
    if(debug==1)fprintf(FILE_PLOT," %i\n", ni_n_b);
    ni_n_b = 301;
    if(debug==1)fprintf(FILE_PLOT," %i\n", ni_n_b);
    for (i_n_b=1; i_n_b<(ni_n_b+1); i_n_b++)
    {
        n_b = 1.e-12*pow(f_n_b,(double)(i_n_b-1));
        if(debug==1)fprintf(FILE_PLOT," %14.7e\n", n_b);
    }
    if(debug==1) fclose(FILE_PLOT);
    /* T */
    if(debug==1)FILE_PLOT = fopen(changestring("t",cnsfile,".dat"),"w");
    ni_t = 1;
    if(debug==1)fprintf(FILE_PLOT," %i\n", ni_t);
    if (ir_t == 0)
    {
        ni_t = 81;
        if(debug==1)fprintf(FILE_PLOT," %i\n", ni_t);
        for (i_t=1; i_t<(ni_t+1); i_t++)
        {
            t = t_ref*pow(f_t,(double)(i_t-1));
            if(debug==1)fprintf(FILE_PLOT," %14.7e\n", t);
        }
    }
    else
    {
        ni_t = 62;
        if(debug==1)fprintf(FILE_PLOT," %i\n", ni_t);
        for (i_t=1; i_t<(ni_t+1); i_t++)
        {
            t = t_ref*sinh(f_t*(double)i_t);
            if(debug==1)fprintf(FILE_PLOT," %14.7e\n", t);
        }
    }
    if(debug==1) fclose(FILE_PLOT);


    return 0;
}
/*****************************************************************************/
int get_cpl(double rho,int icpl)
{
    double x,tmp0,tmp1,tmp2,alpha,beta,f,fp,g,gp,k,p,x0;
    int imf;

    x = rho/RHOREF;

    if (PARA_RMF == 4)
    {
        for (imf=1; imf<5; imf++)
        {
            tmp0 = x+COEFF[imf][4];
            tmp1 = x+COEFF[imf][5];
            tmp2 = 1.+COEFF[imf][3]*QUAD(tmp1);
            CPL[imf][0] = (1.+COEFF[imf][2]*QUAD(tmp0))/tmp2;
            CPL[imf][1] = 2.*(COEFF[imf][2]*tmp0-COEFF[imf][3]*tmp1
                              +COEFF[imf][2]*COEFF[imf][3]*tmp0*tmp1*
                              (COEFF[imf][5]-COEFF[imf][4]))/(RHOREF*QUAD(tmp2));
            tmp2 = COEFF[imf][0]*COEFF[imf][1];
            CPL[imf][0] *= tmp2;
            CPL[imf][1] *= tmp2;
        }
        CPL[5][0] = COEFF[5][0];
        CPL[5][1] = 0.;
    }
    else
    {

        for (imf=1; imf<3; imf++)
        {
            tmp0 = x+COEFF[imf][4];
            tmp1 = QUAD(tmp0);
            tmp2 = 1.+COEFF[imf][3]*tmp1;
            CPL[imf][0] = COEFF[imf][0]*COEFF[imf][1]/tmp2;
            CPL[imf][1] = CPL[imf][0]
                          *2.*(COEFF[imf][2]-COEFF[imf][3])*tmp0/(RHOREF*tmp2);
            CPL[imf][0] *= (1.+COEFF[imf][2]*tmp1);
        }
        CPL[3][0] = COEFF[3][0]*exp(-COEFF[3][1]*(x-1.));
        CPL[3][1] = -COEFF[3][1]*CPL[3][0]/RHOREF;

        CPL[4][0] = COEFF[4][0];
        CPL[4][1] = 0.;

        CPL[5][0] = COEFF[5][0];
        CPL[5][1] = 0.;

        /* high-density modification */
        if (1 == 0)
        {
            x0 = RHOREF;
            x = rho/x0;
            if (x > 1.)
            {
                /* fit of absolute errors */
                /*k = 0.134;
                  p = 0.0314;*/
                /* fit of  relative errors */
                k = 0.159;
                p = 0.0309;
                /*for (imf=1; imf<2; imf++) {*/
                imf = 1;
                {
                    f = CPL[imf][0];
                    fp = CPL[imf][1];
                    if (imf == 1)
                    {
                        alpha = k*(1.-p);
                        beta = k*(1.+p);
                    }
                    else
                    {
                        alpha = k*(1.+p);
                        beta = k*(1.-p);
                    }
                    tmp0 = x-1.;
                    tmp1 = QUAD(tmp0);
                    tmp2 = tmp1*tmp0;
                    tmp0 = 1.+beta*tmp2;
                    g = (1.+alpha*tmp2)/tmp0;
                    gp = 3.*(alpha-beta)*tmp1/(QUAD(tmp0)*x0);
                    CPL[imf][0] = f*g;
                    CPL[imf][1] = f*gp+fp*g;
                }
            }
        }
    }

    if (icpl == 1)
    {
        for (imf=1; imf<N_MF; imf++)
        {
            tmp0         = CPL[imf][0]/MESON[imf].m2;
            CPL2[imf][1] = 2.*tmp0*CPL[imf][1];
            CPL2[imf][0] = tmp0*CPL[imf][0];
        }
    }

    return 0;
}
/*****************************************************************************/
int init_eos (int ipara,int imb,int in_ton, char*cnsfile)
{

    init_parameters();

    init_rmf(ipara);

    /*get_cpl(0.,1);
    if(debug==1)fprintf(myfile,"CPL[1][0] m  %e %e\n", CPL2[1][0], MESON[1].m*HBARC);
    if(debug==1)fprintf(myfile,"CPL[2][0] m  %e %e\n", CPL2[2][0], MESON[2].m*HBARC);
    if(debug==1)fprintf(myfile,"CPL[3][0] m  %e %e\n", CPL2[3][0], MESON[3].m*HBARC);
    if(debug==1)fprintf(myfile,"CPL[4][0] m  %e %e\n", CPL2[4][0], MESON[4].m*HBARC);
    exit(0);*/

    init_nuclei_ame03(in_ton,cnsfile);

    init_ton(in_ton,cnsfile);

    init_particles(imb);

    return 0;
}
/****************************************************************************/
int init_parameters (void)
{

    /* numerical constants */

    PI   = 4.*atan(1.);
    TPI  = 2.*PI;
    FPI  = 4.*PI;
    PI2  = QUAD(PI);
    TPI2 = 2.*PI2;
    RPI  = sqrt(PI);

    /* physical constants */

    HBARC  = 197.3269631;
    HBARC3 = CUBE(HBARC);
    ALPHA  = 0;//1./137.035999679; // Coloumb force switch
    E2     = ALPHA*HBARC;
    G_G    = sqrt(FPI*ALPHA);
    AMU    = 931.494028;
    M_PRO  = 938.272013;
    M_NEU  = 939.565346;

    /*printf(" *** hallo original DD2 masses ***\n");
    M_PRO  = 938.272030;
    M_NEU  = 939.565360;*/

    M_PI0 = 134.9766;
    M_PIP = 139.57018;

    M_K0 = 497.614;
    M_KP = 493.677;

    M_ETA = 547.853;

    M_LAMBDA = 1115.683;
    M_SIGMAP = 1189.37;
    M_SIGMA0 = 1192.642;
    M_SIGMAM = 1197.449;
    M_XI0    = 1314.86;
    M_XIM    = 1321.71;

    U_LAMBDA = -30.;
    U_SIGMA  = -30.;
    U_XI     = -21.;

    R_HYP = 0.83;
    R_PHI = 1.*R_HYP;

    /* conversion factors for units */
    /* nucleon/fm^3 -> g/cm^3 (from eV -> g) */
    CONV[0] = AMU*1.782661758e-36;
    CONV[0] *= 1.e+48;
    /* MeV/fm^3 -> erg/cm^3 (from eV -> J -> erg) */
    CONV[1] = 1.602176487e-19;
    CONV[1] *= 1.e+52;
    /* MeV -> K */
    CONV[2] = 1.1604505e+10;

    return 0;
}
/*****************************************************************************/
int init_rmf (int ipara)
{
    int imf;
    double m_p,m_n,m_omega,m_rho,m_sigma,m_delta,
           rhosat,av,meff,k,j,l,fp,fppfp,
           m_nuc,kf,pf,mu,md,el,rhos,s,v,epsilon,vr,
           c_sigma,c_omega,a,b,c,d,d_old,x,rhosp,mdp,tmp,tmp0,tmp1,tmp2,
           cp_sigma,cpp_sigma,cp_omega,cpp_omega,
           gamma_omega,gamma_sigma,gamma_rho,gamma_delta,
           fp_omega,fpp_omega,fp_sigma,fpp_sigma,c_rho,cp_rho,fp_rho,
           m_phi;

    PARA_RMF = ipara;

    if (ipara == 3)
    {
        NL = 1;
    }
    else
    {
        NL = 0;
    }

    switch (ipara)
    {
    case 2:   /* parametrization DD2 */
    {
        /* original */
        /*m_p   = 938.27203;
          m_n   = 939.56536;*/
        m_p     = M_PRO;
        m_n     = M_NEU;
        m_omega = 783.;
        m_rho   = 763.;
        m_sigma = 546.212459;
        m_delta = 983.;
        /* new TF fit, 2011/10/14 */
        m_sigma = 577.9;
        /* Thomas-Fermi fit 1 ??? */
        /*m_sigma = 586.620187;*/
        /* Thomas-Fermi fit 2 */
        /*m_sigma = 482.612895;*/
        /*m_omega = 559.480364;*/
        rhosat  =   0.149065;
        av      = -16.022477;
        meff    =   0.562526;
        k       = 242.717148;
        j       =  31.670194;
        l       =  55.035224;
        fp      =  -0.141240;
        fppfp   =  -1.066954;
        break;
    }
    case 3:   /* parametrization TM1 */
    {
        m_p     = 938.;
        m_n     = 938.;
        m_omega = 783.;
        m_rho   = 770.;
        m_sigma = 511.19777;
        m_delta = 980.;

        NL_GS = 10.02892;
        NL_GO = 12.61394;
        NL_GR =  4.63219;
        NL_G2 = -7.23247;
        NL_G3 =  0.61833;
        NL_C3 = 71.30747;

        break;
    }
    case 4:   /* parametrization DD-MEdelta*/
    {
        m_p     = M_PRO;
        m_n     = M_NEU;
        m_omega = 783.;
        m_rho   = 763.;
        m_sigma = 566.1577;
        m_delta = 983.;
        break;
    }
    default :   /* parametrization DD */
    {
        m_p     = 939.;
        m_n     = 939.;
        m_omega = 783.;
        m_rho   = 763.;
        m_sigma = 547.204590;
        m_delta = 980.;
        rhosat  =   0.148746;
        av      = -16.021000;
        meff    =   0.564933;
        k       = 239.983000;
        j       =  31.639000;
        l       =  55.976000;
        fp      =  -0.145200;
        fppfp   =  -1.099830;
        break;
    }
    }

    /* phi(1020) */
    m_phi = 1019.455;

    m_nuc = 0.5*(m_p+m_n);

    /* photon */
    MESON[0].m = 0.;
    /* omega meson */
    MESON[1].m = m_omega/HBARC;
    /* sigma meson */
    MESON[2].m = m_sigma/HBARC;
    /* rho meson */
    MESON[3].m = m_rho/HBARC;
    /* delta meson */
    MESON[4].m = m_delta/HBARC;
    /* phi meson */
    MESON[5].m = m_phi/HBARC;

    for (imf=0; imf<N_MF; imf++) MESON[imf].m2 = QUAD(MESON[imf].m);

    if (NL == 1)
    {
        COEFF[1][0] = 1.;
        COEFF[1][1] = 1.;
        COEFF[1][2] = 0.;
        COEFF[1][3] = 0.;
        COEFF[1][4] = 0.;
        COEFF[2][0] = 1.;
        COEFF[2][1] = 1.;
        COEFF[2][2] = 0.;
        COEFF[2][3] = 0.;
        COEFF[2][4] = 0.;
        COEFF[3][0] = 1.;
        COEFF[3][1] = 0.;
        COEFF[4][0] = 0.;
        COEFF[4][1] = 0.;
        COEFF[5][0] = 0.;
        COEFF[5][1] = 0.;
        return 0;
    }

    if (ipara == 4)
    {
        gamma_omega = 12.2904;
        COEFF[1][0] = gamma_omega;
        COEFF[1][1] = 1.4089;
        COEFF[1][2] = 0.1698;
        COEFF[1][3] = 0.3429;
        COEFF[1][4] = 0.9860;
        COEFF[1][5] = 0.9860;
        gamma_sigma = 10.3325;
        COEFF[2][0] = gamma_sigma;
        COEFF[2][1] = 1.3927;
        COEFF[2][2] = 0.1901;
        COEFF[2][3] = 0.3679;
        COEFF[2][4] = 0.9519;
        COEFF[2][5] = 0.9519;
        gamma_rho   =  6.3128;
        COEFF[3][0] = gamma_rho;
        COEFF[3][1] = 1.8877;
        COEFF[3][2] = 0.0651;
        COEFF[3][3] = 0.3469;
        COEFF[3][4] = 0.9417;
        COEFF[3][5] = 0.9737;
        gamma_delta =  7.1520;
        COEFF[4][0] = gamma_delta;
        COEFF[4][1] = 1.5178;
        COEFF[4][2] = 0.3262;
        COEFF[4][3] = 0.6041;
        COEFF[4][4] = 0.4257;
        COEFF[4][5] = 0.5885;
        COEFF[5][0] = 0.;
        COEFF[5][1] = 0.;
        RHOREF = 0.152;
    }
    else
    {

        kf = pow(1.5*PI2*rhosat,0.333333333333);
        pf = HBARC*kf;
        mu = m_nuc+av;
        md = m_nuc*meff;
        el = sqrt(md*md+pf*pf);
        rhos = md*(pf*el-md*md*log((pf+el)/md))/(PI2*HBARC3);
        s = m_nuc-md;
        v = mu-el;
        c_sigma = s/rhos;
        epsilon = mu*rhosat;
        c_omega = 2.*(epsilon-0.75*el*rhosat-0.25*md*rhos
                      -0.5*c_sigma*QUAD(rhos))/QUAD(rhosat);
        vr = v-c_omega*rhosat;
        /*p = 0.25*(el*rhosat-md*rhos)
          -0.5*(c_sigma*QUAD(rhos)-c_omega*QUAD(rhosat))
          +rhosat*vr;*/

        RHOREF = rhosat;
        RHOSREF = rhos;

        x = fppfp;
        d = 0.;
        do
        {
            d_old = d;
            tmp = (-3.-6.*d)/((1.+d)*x);
            tmp = 0.25*(tmp-1.)+0.0625;
            d = sqrt(tmp)-0.25;
        }
        while (fabs(d-d_old) > 1.e-10);

        fp_omega  = fp;
        fpp_omega = fp_omega*x;
        cp_omega = 2.*fp_omega*c_omega/rhosat;
        cpp_omega = 2.*(QUAD(fp_omega)+fpp_omega)*c_omega
                    /QUAD(rhosat);

        cp_sigma = (cp_omega*QUAD(rhosat)-2.*vr)/QUAD(rhos);

        tmp = rhos/md-rhosat/el;
        rhosp = (md/el-3.*cp_sigma*rhos*tmp)/(1.+3.*c_sigma*tmp);

        cpp_sigma = k/9.-QUAD(pf)/(3.*el)-(c_omega+2.*cp_omega*rhosat
                                           +0.5*cpp_omega*QUAD(rhosat))*rhosat;
        cpp_sigma = cpp_sigma+(c_sigma*rhosp*md/el
                               +cp_sigma*rhos*md/el+cp_sigma*rhos*rhosp)*rhosat;
        cpp_sigma = -2.*cpp_sigma/(rhosat*QUAD(rhos));

        fp_sigma  = 0.5*cp_sigma*rhosat/c_sigma;
        fpp_sigma = 0.5*cpp_sigma*QUAD(rhosat)/c_sigma
                    -QUAD(fp_sigma);

        c_rho = 2.*(j-QUAD(pf)/(6.*el))/rhosat;
        mdp = -(cp_sigma*rhos+c_sigma*md/el)
              /(1.+3.*c_sigma*(rhos/md-rhosat/el));
        cp_rho = (l-(1.+QUAD(md)/QUAD(el)
                     -3.*mdp*rhosat*md/QUAD(el))*QUAD(pf)/(6.*el)
                  -1.5*c_rho*rhosat)
                 /(1.5*QUAD(rhosat));

        gamma_omega = MESON[1].m*sqrt(c_omega/HBARC);
        gamma_sigma = MESON[2].m*sqrt(c_sigma/HBARC);
        gamma_rho   = MESON[3].m*sqrt(c_rho/HBARC);
        gamma_delta = 0.;

        fp_rho = 0.5*rhosat*cp_rho/c_rho;

        c = 1./(3.*QUAD(d));
        tmp0 = 1.+d;
        tmp1 = QUAD(tmp0);
        tmp2 = 2.*tmp0/(1.+c*tmp1);
        b = (c*tmp2+fp_omega)/(tmp2-fp_omega*tmp1);
        a = (1.+c*tmp1)/(1.+b*tmp1);
        COEFF[1][0] = gamma_omega;
        COEFF[1][1] = a;
        COEFF[1][2] = b;
        COEFF[1][3] = c;
        COEFF[1][4] = d;

        x = fpp_sigma/fp_sigma;
        if ((x < -3.) || (x >= 0.))
        {
            call_error(1);
        }
        d = 0.;
        do
        {
            d_old = d;
            tmp = (-3.-6.*d)/((1.+d)*x);
            tmp = 0.25*(tmp-1.)+0.0625;
            d = sqrt(tmp)-0.25;
        }
        while (fabs(d-d_old) > 1.e-10);
        x = -2.*(1.+d)/(4.*QUAD(d)+2.*d+1.);
        if ((fp_sigma <= x) || (fp_sigma >= 0.))
        {
            call_error(2);
        }
        c = 1./(3.*QUAD(d));
        tmp0 = 1.+d;
        tmp1 = QUAD(tmp0);
        tmp2 = 2.*tmp0/(1.+c*tmp1);
        b = (c*tmp2+fp_sigma)/(tmp2-fp_sigma*tmp1);
        a = (1.+c*tmp1)/(1.+b*tmp1);
        COEFF[2][0] = gamma_sigma;
        COEFF[2][1] = a;
        COEFF[2][2] = b;
        COEFF[2][3] = c;
        COEFF[2][4] = d;

        COEFF[3][0] = gamma_rho;
        COEFF[3][1] = -fp_rho;

        COEFF[4][0] = gamma_delta;
        COEFF[4][1] = 0.;

        COEFF[5][0] = gamma_omega;
        COEFF[5][1] = 0.;

    }

    get_cpl(0.,0);
    if (CPL[1][0] > 0.)
    {
        LA_OMEGA = MESON[1].m2/CPL[1][0];
    }
    else
    {
        LA_OMEGA = 0.;
    }
    if (CPL[3][0] > 0.)
    {
        LA_RHO   = MESON[3].m2/CPL[3][0];
    }
    else
    {
        LA_RHO = 0.;
    }

    /* no interaction */
    if (ipara == 0)
    {
        if(debug==1)fprintf(myfile," no interaction\n");
        COEFF[1][0] = COEFF[2][0] = COEFF[3][0] = COEFF[4][0] = COEFF[5][0] = 0.;
        LA_OMEGA = LA_RHO = 0.;
    }

    /*
    if(debug==1)fprintf(myfile," %f %f %f %f %f %f\n",
     COEFF[1][0], COEFF[1][1], COEFF[1][2], COEFF[1][3], COEFF[1][4], COEFF[1][5]);
    if(debug==1)fprintf(myfile," %f %f %f %f %f %f\n",
     COEFF[2][0], COEFF[2][1], COEFF[2][2], COEFF[2][3], COEFF[2][4], COEFF[2][5]);
    if(debug==1)fprintf(myfile," %f %f\n", COEFF[3][0], COEFF[4][0]);
    */

    /*printf(" constant couplings\n");
    COEFF[1][2] = COEFF[1][3] = COEFF[1][4] = 0.;
    COEFF[2][2] = COEFF[2][3] = COEFF[2][4] = 0.;
    COEFF[3][1] = COEFF[4][1] = COEFF[5][1] = 0.;*/
    /*COEFF[3][0] = COEFF[4][0] = COEFF[5][0] = 0.;*/
    /*COEFF[1][0] *= 0.5;
      COEFF[2][0] *= 0.5;*/

    return 0;
}
/*****************************************************************************/
int init_ton (int in_ton, char *cnsfile)
{
    int ip,iq,itmp[3],idx;
    float tmp;
    double m_p,m_n;
    FILE *file_bea;

    m_p = M_PRO/HBARC;
    m_n = M_NEU/HBARC;

    for(ip=0; ip<N_AME11; ip++)
    {
        NUC[ip].n   = 0;
        NUC[ip].z   = 0;
        NUC[ip].a   = 0;
        NUC[ip].b   = 0.;
        NUC[ip].bea = 0.;
        NUC[ip].g0  = 0.;
        NUC[ip].g   = 0.;
        NUC[ip].m0  = 0.;
        NUC[ip].m   = 0.;
    }

    if (in_ton == 1)
    {

        file_bea = fopen(changestring("ame11",cnsfile,".in"),"r");

        iq = -1;
        for(ip=0; ip<N_AME11; ip++)
        {
            fscanf(file_bea," %i %i %i %f",
                   &itmp[0], &itmp[1], &itmp[2], &tmp);
            /*printf(" ip n z a %i %i %i %i\n", ip, itmp[0], itmp[1], itmp[2]);*/

            /*tmp = 0.;*/

            /* hallo */
            idx = 0;

            /* remove light clusters */
            if (1 == 1)
            {
                /* no 2H */
                if ((itmp[0] == 1) && (itmp[1] == 1)) idx = 1;
                /* no 3H */
                if ((itmp[0] == 2) && (itmp[1] == 1)) idx = 1;
                /* no 3He */
                if ((itmp[0] == 1) && (itmp[1] == 2)) idx = 1;
                /* no 4He */
                if ((itmp[0] == 2) && (itmp[1] == 2)) idx = 1;
            }

            if (idx == 1)
            {
                iq += 1;
                NUC[iq].n   = itmp[0];
                NUC[iq].z   = itmp[1];
                NUC[iq].a   = itmp[2];
                NUC[iq].b   = (double)(NUC[iq].z-NUC[iq].n)/(double)NUC[iq].a;
                NUC[iq].bea = (double)tmp/HBARC;
                NUC[iq].m0  = (double)NUC[iq].n*m_n
                              +(double)NUC[iq].z*m_p-(double)NUC[iq].a*NUC[iq].bea;

                if ((NUC[iq].a%2) == 1)
                {
                    NUC[iq].g0 = 2.;
                }
                else
                {
                    NUC[iq].g0 = 1.;
                }

                /*printf(" %i %i %i %i %f %f %f\n",
                       iq, NUC[iq].n, NUC[iq].z, NUC[iq].a,
                       NUC[iq].bea*HBARC, NUC[iq].m0*HBARC, NUC[iq].g0);*/
            }
        }
        N_NUC = iq+1;
        /*printf(" N_NUC = %i\n", N_NUC);*/

        fclose(file_bea);

    }

    return 0;
}
/*****************************************************************************/
int init_nuclei_ame03 (int in_ton,char *cnsfile)
{
    int ip,iq,itmp[3],tjp[120][230];
    float tmp;
    double m_p,m_n;
    FILE *file_bea,*file_sp;

    m_p = M_PRO/HBARC;
    m_n = M_NEU/HBARC;

    for(ip=0; ip<N_AME03; ip++)
    {
        NUCLEUS[ip].a   = 0;
        NUCLEUS[ip].n   = 0;
        NUCLEUS[ip].z   = 0;
        NUCLEUS[ip].s   = 0;
        NUCLEUS[ip].st  = 0;
        NUCLEUS[ip].in  = 0;
        NUCLEUS[ip].g   = 0.;
        NUCLEUS[ip].be  = 0.;
        NUCLEUS[ip].m   = 0.;
        NUCLEUS[ip].fs  = 0.;
        NUCLEUS[ip].fo  = 0.;
        NUCLEUS[ip].fr  = 0.;
        NUCLEUS[ip].fd  = 0.;
        NUCLEUS[ip].rms = 0.;
        NUCLEUS[ip].ag  = 0.;
    }

    /*if (in_ton == 1) {*/
    if (0 == 1)
    {

        file_bea = fopen(changestring("ame03",cnsfile,".in"),"r");

        for(ip=0; ip<N_AME03; ip++)
        {
            fscanf(file_bea," %i %i %i %f",
                   &itmp[0], &itmp[1], &itmp[2], &tmp);
            /*printf(" ip a n z %i %i %i %i\n", ip, itmp[2], itmp[0], itmp[1]);*/
            NUCLEUS[ip].a   = itmp[2];
            NUCLEUS[ip].n   = itmp[0];
            NUCLEUS[ip].z   = itmp[1];
            NUCLEUS[ip].st  = NUCLEUS[ip].a%2;
            NUCLEUS[ip].in  = 1;
            NUCLEUS[ip].g   = (double)NUCLEUS[ip].st+1.; /* to be corrected below */
            NUCLEUS[ip].be  = (double)NUCLEUS[ip].a*(double)tmp/HBARC;
            NUCLEUS[ip].m   = (double)NUCLEUS[ip].n*m_n
                              +(double)NUCLEUS[ip].z*m_p-NUCLEUS[ip].be;
            /*
            if (ip < 10)if(debug==1)fprintf(myfile,"ip   a z n   %i   %i %i %i   %f\n", ip,
            	  NUCLEUS[ip].a, NUCLEUS[ip].z, NUCLEUS[ip].n,
            	  NUCLEUS[ip].be*HBARC);
            */
            /*if (ip < 10)if(debug==1)fprintf(myfile,"ip   a z n   %i   %i %i %i   %f\n", ip,
            	  NUCLEUS[ip].a, NUCLEUS[ip].z, NUCLEUS[ip].n,
            	  tmp);*/
        }

        fclose(file_bea);

    }
    else
    {
        /* 1n */
        NUCLEUS[0].a   = 1;
        NUCLEUS[0].n   = 1;
        NUCLEUS[0].z   = 0;
        NUCLEUS[0].be  = 0.;
        /* 1H */
        NUCLEUS[1].a   = 1;
        NUCLEUS[1].n   = 0;
        NUCLEUS[1].z   = 1;
        NUCLEUS[1].be  = 0.;
        /* 2H */
        NUCLEUS[2].a   = 2;
        NUCLEUS[2].n   = 1;
        NUCLEUS[2].z   = 1;
        NUCLEUS[2].be  = (double)NUCLEUS[2].a*1.112283/HBARC;
        /* 3H */
        NUCLEUS[3].a   = 3;
        NUCLEUS[3].n   = 2;
        NUCLEUS[3].z   = 1;
        NUCLEUS[3].be  = (double)NUCLEUS[3].a*2.827266/HBARC;
        /* 4H */
        NUCLEUS[4].a   = 4;
        NUCLEUS[4].n   = 3;
        NUCLEUS[4].z   = 1;
        NUCLEUS[4].be  = (double)NUCLEUS[4].a*1.400351/HBARC;
        /* 5H */
        NUCLEUS[5].a   = 5;
        NUCLEUS[5].n   = 4;
        NUCLEUS[5].z   = 1;
        NUCLEUS[5].be  = (double)NUCLEUS[5].a*1.336360/HBARC;
        /* 6H */
        NUCLEUS[6].a   = 6;
        NUCLEUS[6].n   = 5;
        NUCLEUS[6].z   = 1;
        NUCLEUS[6].be  = (double)NUCLEUS[6].a*0.963633/HBARC;
        /* 7H */
        NUCLEUS[7].a   = 7;
        NUCLEUS[7].n   = 6;
        NUCLEUS[7].z   = 1;
        NUCLEUS[7].be  = (double)NUCLEUS[7].a*0.940000/HBARC;
        /* 3He */
        NUCLEUS[8].a   = 3;
        NUCLEUS[8].n   = 1;
        NUCLEUS[8].z   = 2;
        NUCLEUS[8].be  = (double)NUCLEUS[8].a*2.572681/HBARC;
        /* 4He */
        NUCLEUS[9].a   = 4;
        NUCLEUS[9].n   = 2;
        NUCLEUS[9].z   = 2;
        NUCLEUS[9].be  = (double)NUCLEUS[9].a*7.073915/HBARC;

        for(ip=0; ip<10; ip++)
        {
            NUCLEUS[ip].st  = NUCLEUS[ip].a%2;
            NUCLEUS[ip].g   = (double)NUCLEUS[ip].st+1.; /* to be corrected below */
            NUCLEUS[ip].m   = (double)NUCLEUS[ip].n*m_n
                              +(double)NUCLEUS[ip].z*m_p-NUCLEUS[ip].be;
            /*if (ip < 10)if(debug==1)fprintf(myfile,"ip   a z n   %i   %i %i %i   %f\n", ip,
            	  NUCLEUS[ip].a, NUCLEUS[ip].z, NUCLEUS[ip].n,
            	  NUCLEUS[ip].be*HBARC);*/
        }

    }

    for (ip=0; ip<120; ip++)
    {
        for (iq=0; iq<230; iq++) tjp[ip][iq] = 0;
    }

    if (in_ton == 1)
    {

        file_sp = fopen(changestring("spin",cnsfile,".in"),"r");

        for(ip=0; ip<N_SPIN; ip++)
        {
            fscanf(file_sp," %i %i %i",
                   &itmp[0], &itmp[1], &itmp[2]);
            if ((itmp[0] < 120) && (itmp[1] < 230))
            {
                tjp[itmp[0]][itmp[1]] = itmp[2]+1;
            }
            /*
            if (ip < 16) {
            if(debug==1)fprintf(myfile," %i %i %i\n", itmp[0], itmp[1], tjp[itmp[0]][itmp[1]]);
            }
            */
        }

        fclose(file_sp);

    }
    else
    {
        tjp[0][1] = 2;
        tjp[1][0] = 2;
        tjp[1][1] = 3;
        tjp[1][2] = 2;
        tjp[2][1] = 2;
        tjp[2][2] = 1;
        tjp[2][4] = 1;
    }

    for (ip=0; ip<N_AME03; ip++)
    {
        if ((NUCLEUS[ip].n%2 == 1) || (NUCLEUS[ip].z%2 == 1))
        {
            if (tjp[NUCLEUS[ip].z][NUCLEUS[ip].n] > 0)
            {
                NUCLEUS[ip].g = (double)tjp[NUCLEUS[ip].z][NUCLEUS[ip].n];
            }
        }
    }

    return 0;
}
/*****************************************************************************/
int init_particles(int imb)
{
    int ip,icl;
    double tmp,a,b,c,d,g,l,m,b2,g2,tmp1,tmp2,g_p,g_m,fo_lambda,fo_sigma,fo_xi,
           fs_lambda,fs_sigma,fs_xi;

    tmp = sqrt(2./3.);

    for (ip=0; ip<N_PART; ip++)
    {
        PARTICLE[ip].a  = 0;
        PARTICLE[ip].n  = 0;
        PARTICLE[ip].z  = 0;
        PARTICLE[ip].s  = 0;
        PARTICLE[ip].st = 0;
        PARTICLE[ip].g  = 0.;
        PARTICLE[ip].g0 = 0.;
        PARTICLE[ip].be = 0.;
        PARTICLE[ip].m  = 0.;
        PARTICLE[ip].fs = 0.;
        PARTICLE[ip].fo = 0.;
        PARTICLE[ip].fr = 0.;
        PARTICLE[ip].fd = 0.;
        PARTICLE[ip].fp = 0.;
        PARTICLE[ip].rms= 0.;
        PARTICLE[ip].ag = 0.;
        PARTICLE[ip].mu = 0.;
    }

    /* particle statistics:
       -2: Bose-Einstein, massless
       -1: Bose-Einstein, boson identical with antiboson
       0: Bose-Einstein, boson not identical with antiboson
       1: Fermi-Dirac
       2: Maxwell-Boltzmann, non-relativistic
       3: Maxwell-Boltzmann, relativistic
    */

    /* proton */
    ip = 0;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_PRO/HBARC;
    PARTICLE[ip].fs =  1.;
    PARTICLE[ip].fo =  1.;
    PARTICLE[ip].fr =  1.;
    PARTICLE[ip].fd =  1.;
    PARTICLE[ip].rms=  0.; /*0.64;*/
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;
    NUCLEUS[1].in = 0;

    /* neutron */
    ip = 1;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_NEU/HBARC;
    PARTICLE[ip].fs =  1.;
    PARTICLE[ip].fo =  1.;
    PARTICLE[ip].fr = -1.;
    PARTICLE[ip].fd = -1.;
    PARTICLE[ip].rms=  0.; /*0.64;*/
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;
    NUCLEUS[0].in = 0;

    /*PARTICLE[0].st =  2;
      PARTICLE[1].st =  2;*/

    /*PARTICLE[0].g = 0.000001;
      PARTICLE[1].g = 0.000001;*/

    /* electron */
    ip = 2;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  0.510998910/HBARC;

    /* muon */
    ip = 3;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  105.6583668/HBARC;

    /* tauon */
    ip = 4;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  1776.82/HBARC;

    /* electron neutrino */
    ip = 5;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = 1.;

    /* muon neutrino */
    ip = 6;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = 1.;

    /* tau neutrino */
    ip = 7;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = 1.;

    /* photon */
    ip = 8;
    PARTICLE[ip].st = 4;
    PARTICLE[ip].g  = 2.;

    /* deuteron */
    ip = 9;
    PARTICLE[ip].a  = 2;
    PARTICLE[ip].n  = 1;
    PARTICLE[ip].z  = 1;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = 3.;
    PARTICLE[ip].be = (double)PARTICLE[ip].a*1.112283/HBARC;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m-PARTICLE[ip].be;
    PARTICLE[ip].fs = 2.;
    PARTICLE[ip].fo = 2.;
    PARTICLE[ip].fr = 0.;
    PARTICLE[ip].fd = 0.;
    PARTICLE[ip].rms= 1.96;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;
    NUCLEUS[2].in = 0;

    /* triton */
    ip = 10;
    PARTICLE[ip].a  = 3;
    PARTICLE[ip].n  = 2;
    PARTICLE[ip].z  = 1;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = 2.;
    PARTICLE[ip].be = (double)PARTICLE[ip].a*2.827266/HBARC;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m-PARTICLE[ip].be;
    PARTICLE[ip].m  = NUCLEUS[3].m;
    PARTICLE[ip].fs = 3.;
    PARTICLE[ip].fo = 3.;
    PARTICLE[ip].fr = -1.;
    PARTICLE[ip].fd = -1.;
    PARTICLE[ip].rms= 1.59;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;
    NUCLEUS[3].in = 0;

    /* helion */
    ip = 11;
    PARTICLE[ip].a  = 3;
    PARTICLE[ip].n  = 1;
    PARTICLE[ip].z  = 2;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = 2.;
    PARTICLE[ip].be = (double)PARTICLE[ip].a*2.572681/HBARC;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m-PARTICLE[ip].be;
    PARTICLE[ip].m  = NUCLEUS[8].m;
    PARTICLE[ip].fs = 3.;
    PARTICLE[ip].fo = 3.;
    PARTICLE[ip].fr = 1.;
    PARTICLE[ip].fd = 1.;
    PARTICLE[ip].rms= 1.76;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;
    NUCLEUS[8].in = 0;

    /* alpha */
    ip = 12;
    PARTICLE[ip].a  = 4;
    PARTICLE[ip].n  = 2;
    PARTICLE[ip].z  = 2;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = 1.;
    PARTICLE[ip].be = (double)PARTICLE[ip].a*7.073915/HBARC;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m-PARTICLE[ip].be;
    PARTICLE[ip].m  = NUCLEUS[9].m;
    PARTICLE[ip].fs = 4.;
    PARTICLE[ip].fo = 4.;
    PARTICLE[ip].fr = 0.;
    PARTICLE[ip].fd = 0.;
    PARTICLE[ip].rms= 1.45;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;
    NUCLEUS[9].in = 0;

    /* np 3S1 */
    ip = 13;
    PARTICLE[ip].a  = 2;
    PARTICLE[ip].n  = 1;
    PARTICLE[ip].z  = 1;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = -3.;
    PARTICLE[ip].g0 = 3.;
    PARTICLE[ip].be = PARTICLE[9].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 2.;
    PARTICLE[ip].fo = 2.;
    PARTICLE[ip].fr = 0.;
    PARTICLE[ip].fd = 0.;
    PARTICLE[ip].rms= PARTICLE[9].rms;
    PARTICLE[ip].ag = PARTICLE[9].ag;

    /* np 1S0 */
    ip = 14;
    PARTICLE[ip].a  = 2;
    PARTICLE[ip].n  = 1;
    PARTICLE[ip].z  = 1;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = 1.;
    PARTICLE[ip].g0 = 1.;
    PARTICLE[ip].be = PARTICLE[9].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 2.;
    PARTICLE[ip].fo = 2.;
    PARTICLE[ip].fr = 0.;
    PARTICLE[ip].fd = 0.;
    PARTICLE[ip].rms= PARTICLE[9].rms;
    PARTICLE[ip].ag = PARTICLE[9].ag;

    /* nn 1S0 */
    ip = 15;
    PARTICLE[ip].a  = 2;
    PARTICLE[ip].n  = 2;
    PARTICLE[ip].z  = 0;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = 1.;
    PARTICLE[ip].g0 = 1.;
    PARTICLE[ip].be = PARTICLE[9].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 2.;
    PARTICLE[ip].fo = 2.;
    PARTICLE[ip].fr = -2.;
    PARTICLE[ip].fd = -2.;
    /*
    PARTICLE[ip].fs = -2.;
    PARTICLE[ip].fo = -2.;
    PARTICLE[ip].fr = 2.;
    PARTICLE[ip].fd = 2.;
    */
    PARTICLE[ip].rms= PARTICLE[9].rms;
    PARTICLE[ip].ag = PARTICLE[9].ag;

    /* pp 1S0 */
    ip = 16;
    PARTICLE[ip].a  = 2;
    PARTICLE[ip].n  = 0;
    PARTICLE[ip].z  = 2;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = 1.;
    PARTICLE[ip].g0 = 1.;
    PARTICLE[ip].be = PARTICLE[9].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 2.;
    PARTICLE[ip].fo = 2.;
    PARTICLE[ip].fr = 2.;
    PARTICLE[ip].fd = 2.;
    PARTICLE[ip].rms= PARTICLE[9].rms;
    PARTICLE[ip].ag = PARTICLE[9].ag;

    /* triton continuum */
    ip = 17;
    PARTICLE[ip].a  = 3;
    PARTICLE[ip].n  = 2;
    PARTICLE[ip].z  = 1;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = -2.;
    PARTICLE[ip].be = PARTICLE[10].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 3.;
    PARTICLE[ip].fo = 3.;
    PARTICLE[ip].fr = -1.;
    PARTICLE[ip].fd = -1.;
    PARTICLE[ip].rms= 1.59;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;

    /* helion continuum */
    ip = 18;
    PARTICLE[ip].a  = 3;
    PARTICLE[ip].n  = 1;
    PARTICLE[ip].z  = 2;
    PARTICLE[ip].st = 1;
    PARTICLE[ip].g  = -2.;
    PARTICLE[ip].be = PARTICLE[11].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 3.;
    PARTICLE[ip].fo = 3.;
    PARTICLE[ip].fr = 1.;
    PARTICLE[ip].fd = 1.;
    PARTICLE[ip].rms= 1.76;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;

    /* alpha continuum */
    ip = 19;
    PARTICLE[ip].a  = 4;
    PARTICLE[ip].n  = 2;
    PARTICLE[ip].z  = 2;
    PARTICLE[ip].st = 0;
    PARTICLE[ip].g  = -1.;
    PARTICLE[ip].be = PARTICLE[12].be;
    PARTICLE[ip].m  = (double)PARTICLE[ip].z*PARTICLE[0].m
                      +(double)PARTICLE[ip].n*PARTICLE[1].m;
    PARTICLE[ip].fs = 4.;
    PARTICLE[ip].fo = 4.;
    PARTICLE[ip].fr = 0.;
    PARTICLE[ip].fd = 0.;
    PARTICLE[ip].rms= 1.45;
    PARTICLE[ip].ag = PARTICLE[ip].rms*tmp;

    /* hyperons */
    tmp1 = RHOREF;
    get_cpl(tmp1,1);
    /*printf(" RHO RHOS %f %f\n", RHOREF, RHOSREF);
     if(debug==1)fprintf(myfile," CPL2 %f %f\n", CPL2[1][0], CPL2[2][0]);*/
    fo_lambda = 2.*R_HYP/3.;
    fo_sigma  = 2.*R_HYP/3.;
    fo_xi     = R_HYP/3.;
    fs_lambda = (fo_lambda*(CPL2[1][0]+0.5*CPL2[1][1]*RHOREF)*RHOREF
                 -U_LAMBDA/HBARC)/((CPL2[2][0]+0.5*CPL2[2][1]*RHOSREF)*RHOSREF);
    fs_sigma  = (fo_sigma*(CPL2[1][0]+0.5*CPL2[1][1]*RHOREF)*RHOREF
                 -U_SIGMA/HBARC)/((CPL2[2][0]+0.5*CPL2[2][1]*RHOSREF)*RHOSREF);
    fs_xi     = (fo_xi*(CPL2[1][0]+0.5*CPL2[1][1]*RHOREF)*RHOREF
                 -U_XI/HBARC)/((CPL2[2][0]+0.5*CPL2[2][1]*RHOSREF)*RHOSREF);
    /*
    if(debug==1)fprintf(myfile," fo_lambda %f\n", fo_lambda);
    if(debug==1)fprintf(myfile," fo_sigma  %f\n", fo_sigma);
    if(debug==1)fprintf(myfile," fo_xi     %f\n", fo_xi);
    if(debug==1)fprintf(myfile," fs_lambda %f\n", fs_lambda);
    if(debug==1)fprintf(myfile," fs_sigma  %f\n", fs_sigma);
    if(debug==1)fprintf(myfile," fs_xi     %f\n", fs_xi);
    exit(0);
    */

    /* Lambda */
    ip = 20;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  0;
    PARTICLE[ip].s  = -1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_LAMBDA/HBARC;
    PARTICLE[ip].fs =  fs_lambda;
    PARTICLE[ip].fo =  fo_lambda;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp = -sqrt(2.)*R_PHI/3.;

    /* Sigma+ */
    ip = 21;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].s  = -1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_SIGMAP/HBARC;
    PARTICLE[ip].fs =  fs_sigma;
    PARTICLE[ip].fo =  fo_sigma;
    PARTICLE[ip].fr =  2.*R_HYP;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp = -sqrt(2.)*R_PHI/3.;

    /* Sigma0 */
    ip = 22;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  0;
    PARTICLE[ip].s  = -1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_SIGMA0/HBARC;
    PARTICLE[ip].fs =  fs_sigma;
    PARTICLE[ip].fo =  fo_sigma;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp = -sqrt(2.)*R_PHI/3.;

    /* Sigma- */
    ip = 23;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  = -1;
    PARTICLE[ip].s  = -1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_SIGMAM/HBARC;
    PARTICLE[ip].fs =  fs_sigma;
    PARTICLE[ip].fo =  fo_sigma;
    PARTICLE[ip].fr = -2.*R_HYP;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp = -sqrt(2.)*R_PHI/3.;

    /* Xi0 */
    ip = 24;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  0;
    PARTICLE[ip].s  = -1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_XI0/HBARC;
    PARTICLE[ip].fs =  fs_xi;
    PARTICLE[ip].fo =  fo_xi;
    PARTICLE[ip].fr =  R_HYP;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp = -sqrt(2.)*R_PHI/3.;

    /* Xi- */
    ip = 25;
    PARTICLE[ip].a  =  1;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  = -1;
    PARTICLE[ip].s  = -1;
    PARTICLE[ip].st =  1;
    PARTICLE[ip].g  =  2.;
    PARTICLE[ip].m  =  M_XIM/HBARC;
    PARTICLE[ip].fs =  fs_xi;
    PARTICLE[ip].fo =  fo_xi;
    PARTICLE[ip].fr = -R_HYP;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp = -sqrt(2.)*R_PHI/3.;

    /* thermal mesons: presently rel. Maxwell-Boltzmann statistics */

    /* pi0 */
    ip = 26;
    PARTICLE[ip].a  =  0;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  0;
    PARTICLE[ip].s  =  0;
    PARTICLE[ip].st =  3;
    PARTICLE[ip].g  =  0.5; /* p = anti-p */
    PARTICLE[ip].m  =  M_PI0/HBARC;
    PARTICLE[ip].fs =  0.;
    PARTICLE[ip].fo =  0.;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp =  0.;

    /* pi+(-) */
    ip = 27;
    PARTICLE[ip].a  =  0;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].s  =  0;
    PARTICLE[ip].st =  3;
    PARTICLE[ip].g  =  1.;
    PARTICLE[ip].m  =  M_PIP/HBARC;
    PARTICLE[ip].fs =  0.;
    PARTICLE[ip].fo =  0.;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp =  0.;

    /* K0(K0bar) */
    ip = 28;
    ip = 28;
    PARTICLE[ip].a  =  0;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  0;
    PARTICLE[ip].s  =  1;
    PARTICLE[ip].st =  3;
    PARTICLE[ip].g  =  1.;
    PARTICLE[ip].m  =  M_K0/HBARC;
    PARTICLE[ip].fs =  0.;
    PARTICLE[ip].fo =  0.;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp =  0.;

    /* K+(-) */
    ip = 29;
    PARTICLE[ip].a  =  0;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  1;
    PARTICLE[ip].s  =  1;
    PARTICLE[ip].st =  3;
    PARTICLE[ip].g  =  1.;
    PARTICLE[ip].m  =  M_KP/HBARC;
    PARTICLE[ip].fs =  0.;
    PARTICLE[ip].fo =  0.;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp =  0.;

    /* eta */
    ip = 30;
    PARTICLE[ip].a  =  0;
    PARTICLE[ip].n  =  0;
    PARTICLE[ip].z  =  0;
    PARTICLE[ip].s  =  0;
    PARTICLE[ip].st =  3;
    PARTICLE[ip].g  =  0.5; /* p = anti-p */
    PARTICLE[ip].m  =  M_ETA/HBARC;
    PARTICLE[ip].fs =  0.;
    PARTICLE[ip].fo =  0.;
    PARTICLE[ip].fr =  0.;
    PARTICLE[ip].fd =  0.;
    PARTICLE[ip].fp =  0.;

    /* clusters */
    /* use Maxwell-Boltzmann statistics for clusters */
    if (imb == 1)
    {
        for (ip=9; ip<20; ip++)
            PARTICLE[ip].st = 2;
    }
    if (imb == 2)
    {
        for (ip=9; ip<20; ip++)
            PARTICLE[ip].st = 3;
    }

    m = 0.5*(939.565+938.783)/HBARC;
    /* Pauli shifts of binding energies in Gaussian/Jastrow approach */
    for (icl=0; icl<4; icl++)
    {
        switch (icl)
        {
        case 0:   /* 2H */
        {
            a = 2.;
            b = 0.625;
            c = 2.;
            d = 0.;
            g = 0.753;
            l = -3677.2/HBARC;
            break;
        }
        case 1:   /* 3H */
        {
            a = 3.;
            b = 0.889;
            c = 24.;
            d = 14.;
            g = 1.083;
            l = -1670.0/HBARC;
            break;
        }
        case 2:   /* 3He */
        {
            a = 3.;
            b = 0.804;
            c = 24.;
            d = 14.;
            g = 0.960;
            l = -1957.5/HBARC;
            break;
        }
        case 3:   /* 4He */
        {
            a = 4.;
            b = 1.034;
            c = 16.;
            d = 10.;
            g = 1.152;
            l = -1449.6/HBARC;
            break;
        }
        default:
        {
            a = b = c = d = l = 0.;
            g = 1.;
        }
        }
        A_SC[icl] = 0.;
        R_SC[icl] = 0.;
        CL_A[icl] = 0.;
        CL_B[icl] = 0.;
        CL_C[icl] = 0.;
        CL_D[icl] = 0.;
        CL_MU[icl] = 0.;
        b2 = QUAD(b);
        g2 = QUAD(g);
        tmp = b2+2.*g2;
        DEPG[0][icl] = -a*(a-1.)*l*b*b2*CUBE(g2)/(sqrt(2.*CUBE(m))*CUBE(tmp));
        DEPG[1][icl] = b2*(b2+d*g2)/(c*m*tmp);
    }
    /* Pauli shifts of binding energies in Jastrow approach */
    /* only for deuteron */
    a = 1.474;
    b = 0.2317;
    DEPJ[0] = 0.25*QUAD(a)/m;
    DEPJ[1] = sqrt(2.)*b/a;
    b = DEPJ[1];
    b2 = QUAD(b);
    tmp = 0.5*sqrt(PI)*exp(b2)*erfc(b)*(1.+2.*b2)/b-1.;
    DEPJ[2] = sqrt(CUBE(TPI/m))*DEPJ[0]/tmp;

    /* shift for A>2 cluster continuum contributions */
    DEPG[0][8]  = DEPG[0][1];
    DEPG[0][9]  = DEPG[0][2];
    DEPG[0][10] = DEPG[0][3];
    DEPG[1][8]  = DEPG[1][1];
    DEPG[1][9]  = DEPG[1][2];
    DEPG[1][10] = DEPG[1][3];

    /* two-body continuum correlations */
    if (1 == 1)
    {
        /*
           scattering lengths and effective ranges from
           R.B. Wiringa, V.G.J. Stoks, R. Schiavilla
           Phys. Rev. C 51, 38 (2001)
           values for v_18 potential without electromagnetic interaction
        */
        /* np 3S1 */
        A_SC[4] = 5.402;
        R_SC[4] = 1.752;
        CL_MU[4] = PARTICLE[0].m*PARTICLE[1].m/(PARTICLE[0].m+PARTICLE[1].m);
        /* np 1S0 */
        A_SC[5] = -23.084;
        R_SC[5] = 2.703;
        CL_MU[5] = CL_MU[4];
        /* nn 1S0 */
        A_SC[6] = -18.818;
        R_SC[6] = 2.834;
        /*A_SC[6] = -23.084;
          R_SC[6] = 2.703;*/
        CL_MU[6] = 0.5*PARTICLE[1].m;
        /* pp 1S0 */
        A_SC[7] = -17.164;
        R_SC[7] = 2.865;
        CL_MU[7] = 0.5*PARTICLE[0].m;
    }
    for (icl=4; icl<8; icl++)
    {
        tmp1 = R_SC[icl]/A_SC[icl];
        if (tmp1 > 0.5)
        {
            if(debug==1)fprintf(myfile," error in init_clusters\n");
            //exit(0);
            return 0;
        }
        tmp2 = sqrt(1.-2.*tmp1);
        CL_C[icl] = 1.+tmp2;
        CL_D[icl] = 1.-tmp2;
        CL_A[icl] = CL_C[icl]-tmp1;
        CL_B[icl] = CL_D[icl]-tmp1;
        tmp1 = 0.5*QUAD(A_SC[icl]);
        CL_A[icl] *= tmp1;
        CL_B[icl] *= tmp1;
        tmp1 = -0.5*A_SC[icl];
        CL_C[icl] *= tmp1;
        CL_D[icl] *= tmp1;
        /*printf(" %i %f %f %f %f\n",
        icl, CL_A[icl], CL_B[icl], CL_C[icl], CL_D[icl]);*/
    }

    if (1 == 0)   /* shifted to init_be_cl */
    {
        get_cpl(0.,0);
        g_p = g_m = QUAD(CPL[2][0])/MESON[2].m2-QUAD(CPL[1][0])/MESON[1].m2;
        tmp1 = QUAD(CPL[4][0])/MESON[4].m2-QUAD(CPL[3][0])/MESON[3].m2;
        /*tmp1 *= -1.;*/
        g_p += tmp1;
        g_m -= tmp1;
        /* pp 1S0 */
        tmp1 = sqrt(PARTICLE[0].m/TPI);
        GGMM[7] = g_p*CUBE(tmp1)/sqrt(2.);
        /* nn 1S0 */
        tmp1 = sqrt(PARTICLE[1].m/TPI);
        GGMM[6] = g_p*CUBE(tmp1)/sqrt(2.);
        tmp1 = sqrt(2.*PARTICLE[0].m*PARTICLE[1].m/
                    (TPI*(PARTICLE[0].m+PARTICLE[1].m)));
        tmp2 = CUBE(tmp1)/sqrt(2.);
        /* pn 1S0 */
        GGMM[5] = g_p*tmp2;
        /* pn 3S1 */
        GGMM[4] = (2.*g_m-g_p)*tmp2;
        /*printf(" G+ = %e\n", g_p);
         if(debug==1)fprintf(myfile," G- = %e\n", g_m);*/
    }

    /*
     if(debug==1)fprintf(myfile," PP1S0 %e %e\n", GGMM[7], HBARC/QUAD(GGMM[7]));
     if(debug==1)fprintf(myfile," NN1S0 %e %e\n", GGMM[6], HBARC/QUAD(GGMM[6]));
     if(debug==1)fprintf(myfile," PN1S0 %e %e\n", GGMM[5], HBARC/QUAD(GGMM[5]));
     if(debug==1)fprintf(myfile," PN3S1 %e %e\n", GGMM[4], HBARC/QUAD(GGMM[4]));
    */

    return 0;
}
/*****************************************************************************/
double experfc(double z)
{
    double z2,z4,z6,res;

    if (z > 0.)
    {
        z2 = QUAD(z);
        if (z < 15.)
        {
            res = exp(z2)*erfc(z);
        }
        else
        {
            z4 = QUAD(z2);
            z6 = z2*z4;
            res = (1.-0.5/z2+0.75/z4-1.875/z6)/(sqrt(PI)*z);
        }
    }
    else
    {
        if(debug==1)fprintf(myfile," error in subroutine experfc ");
        // exit(0);
        return 0;
    }

    return res;
}
/*****************************************************************************/
void call_error (int err)
{

    if(debug==1)fprintf(myfile,"\n   *** ERROR %i ***\n", err);
    if(debug==1)fprintf(myfile,"   *** PROGRAM TERMINATED *** \n\n");

// exit(0);
}
/*****************************************************************************/
int get_nuc_table_old(int inse,int in_g,int in_e,int in_mu,int in_tau,int in_nu,
                      int *in_cl,int in_hyp,int in_tmes,char *cnsfile)
{
    int iwr,itmp,icl,iph,in_cell,inuc,nn_b,in_b,ic,in_ton;
    float tmp;
    double t,n_b,y_q,n_b_min,n_b_max,dn_b,a,z;

    iwr = 0;

    IN_R = 0;
    inse = 0;
    in_g = in_mu = in_tau = in_nu = in_hyp = in_ton = 0;
    for (icl=0; icl<N_CL; icl++) in_cl[icl] = 0;
    IN_PH[0] = IN_PH[1] = 1;
    for (iph=2; iph<N_PH; iph++) IN_PH[iph] = 0;
    in_cell = 1;
    inuc = 1;

    /* read data file nuc_table.in */
    FILE_IN = fopen(changestring("nuc_table",cnsfile,".in"),"r");
    fscanf(FILE_IN,"%i", &itmp);
    a = (double)itmp;
    fscanf(FILE_IN,"%i", &itmp);
    z = (double)itmp;
    fscanf(FILE_IN,"%f", &tmp);
    t = (double)tmp;
    fscanf(FILE_IN,"%f", &tmp);
    n_b_min = (double)tmp;
    fscanf(FILE_IN,"%f", &tmp);
    n_b_max = (double)tmp;
    fscanf(FILE_IN,"%i", &itmp);
    nn_b = itmp;
    fclose(FILE_IN);

    y_q = (double)z/(double)a;

    FILE_OUT_EXTRA  = fopen(changestring("nuc",cnsfile,".dat"),"w");
    if(debug==1)FILE_PLOT       = fopen(changestring("plot",cnsfile,".dat"),"w");
    if(debug==1) FILE_PLOT2      = fopen(changestring("plot2",cnsfile,".dat"),"w");

    /* non-existence of radius solutions before first run */
    for (iph=0; iph<N_PH; iph++)
    {
        IEXR[iph] = 0;
        IMFS[iph] = 0;
    }

    ic = 1;
    if (nn_b > 1)
    {
        dn_b = (n_b_max-n_b_min)/(double)(nn_b-1);
    }
    else
    {
        dn_b = 0.;
    }
    for (in_b=0; in_b<nn_b; in_b++)
    {
        n_b = n_b_min+dn_b*(double)in_b;

        if (n_b > 0.)
        {
            in_e = 1;
        }
        else
        {
            in_e = 0;
        }

        if (ic == 1)
        {
            /***************************************************************/
            get_eos_local(t,n_b,y_q,inse,
                          in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                          in_cell,in_ton,inuc,a,z,cnsfile);
            /***************************************************************/
            if (fabs(AA_NUC-A0_WS) < 1.e-03)
            {
                fprintf(FILE_OUT_EXTRA," %e %e %f %f %e %e %e %e\n",
                        t, n_b, a, z, A0_WS, Z0_WS, BEA_NUC*HBARC, R_WS);
            }
            else
            {
                ic = 0;
            }
        }
    }

    if(debug==1) fclose(FILE_PLOT2);
    if(debug==1) fclose(FILE_PLOT);
    fclose(FILE_OUT_EXTRA);

    return 0;
}
/*****************************************************************************/
double get_nucleus_old(double a,double z)
{
    int iwr,iph;
    double r,tmp1,tmp2,n;

    iwr = 1;

    n = a-z;
    AA_NUC = a;

    init_sol();
    /* homogeneous matter */
    iph = 0;
    if (iwr == 1)if(debug==1)fprintf(myfile," homogeneous matter\n");
    IMFS[iph] = 0;
    iexs_sol(iph);
    r = 0.;
    XCL = 0;
    if (N_B > 0.)
    {
        solve_rmf(iph,r);
    }
    else
    {
        ISOL = 0;
        ICV_RMF = 1;
        NR = 1;
        MU_BB = PARTICLE[1].m;
        MU_QQ = PARTICLE[0].m-PARTICLE[1].m;
        MU_LL = PARTICLE[2].m;
        R_WS = r;
        /*printf(" m_p m_n %e %e\n", PARTICLE[0].m*HBARC, PARTICLE[1].m*HBARC);*/
        F_DENS = z*PARTICLE[0].m+n*PARTICLE[1].m;
    }
    save_sol(iph);
    if (N_B > 0.)
    {
        tmp1 = F_DENS/N_B;
    }
    else
    {
        tmp1 = F_DENS/(double)a;
    }
    if (iwr == 1)if(debug==1)fprintf(myfile," iph F/A %i %f\n", iph, HBARC*tmp1);
    /* nucleus */
    iph = 1;
    if (iwr == 1)if(debug==1)fprintf(myfile," nucleus");
    IMFS[iph] = 0;
    iexs_sol(iph);
    if (N_B > 0.)
    {
        r = pow(3.*(double)a/(FPI*N_B),0.333333333333);
    }
    else
    {
        r = 20.;
    }
    if (iwr == 1)if(debug==1)fprintf(myfile," with cell radius %e fm\n", r);
    XCL = 0;
    if (N_B > 0.)
    {
        fit_radius_nuc(iph,r);
    }
    else
    {
        solve_rmf(iph,r);
    }
    save_sol(iph);
    if (N_B > 0.)
    {
        tmp2 = F_DENS/N_B;
    }
    else
    {
        tmp2 = F_DENS*V_WS/(double)a;
    }
    if (iwr == 1)if(debug==1)fprintf(myfile," iph F/A %i %f\n", iph, HBARC*tmp2);

    if (iwr == 1)if(debug==1)fprintf(myfile," A0 Z0 N0 BE %f %f %f %f\n", A0_WS, Z0_WS, N0_WS, HBARC*(tmp1-tmp2));

    /*
    for (ir=0; ir<NR; ir++) {
      if(debug==1)fprintf(FILE_PLOT," %e %e %e %e\n", RVEC[1][ir],
        DENS[0][ir][2], DENS[1][ir][2], (DENS[0][ir][2]+DENS[1][ir][2]));
    }
    */

    return (tmp1-tmp2);
}
/*****************************************************************************/



/*****************************************************************************/




//MAIN FUNCTION


main(int argc , char * argv[])
{

    /* deifining variables */

    int mpitask_id,number_of_processors;

    int NMZ_MIN , NMZ_MAX , NB_MIN , NB_MAX , T_MIN , T_MAX, A_0_MIN,A_0_MAX;
    int NMZ_SIZE , NB_SIZE , T_SIZE, A_0_SIZE; //range of values of NMZ , NB, T

    int NMZ , T , NB, A_0, i=0, job[4]; //values to be passed to the eos function known as myfunction()

    T_MIN=atoi(argv[1]);
    T_MAX=atoi(argv[2]);
    NB_MIN=atoi(argv[3]);
    NB_MAX=atoi(argv[4]);
    NMZ_MIN=atoi(argv[5]);
    NMZ_MAX=atoi(argv[6]);
    A_0_MIN=atoi(argv[7]);
    A_0_MAX=atoi(argv[8]);
    debug=atoi(argv[9]);

    NB_SIZE=(NB_MAX-NB_MIN)+1;
    NMZ_SIZE=(NMZ_MAX-NMZ_MIN)+1;
    T_SIZE=(T_MAX-T_MIN)+1;
    A_0_SIZE=(A_0_MAX-A_0_MIN)+1;


    long int MAX_SIZE = MAX_SIZE = (NMZ_SIZE)*(NB_SIZE)*(T_SIZE)*(A_0_SIZE);
    /* MPI INITIALIZED */

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &mpitask_id);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);

    if(mpitask_id == 0)
    {
        manager(number_of_processors-1, argc, argv );
    }
    /*  else if(mpitask_id == 1)
      {
          writer();
      }*/
    else
    {
        i = mpitask_id - 1;
        NMZ=(i%NMZ_SIZE)+NMZ_MIN;
        NB=((i/NMZ_SIZE)%NB_SIZE)+NB_MIN;
        T=((i/NB_SIZE)/NMZ_SIZE)%T_SIZE + T_MIN;
        A_0=(((i/NB_SIZE)/NMZ_SIZE)/T_SIZE) + A_0_MIN;


        job[0] = NMZ;
        job[1]= NB;
        job[2] = T;

        job[3] = A_0;

        worker(mpitask_id,MAX_SIZE,job);
    }
    MPI_Finalize();
}

/*function to convert filenames*/

char* changestring(char *pp,char *xx,char*zz) //zz is the filetype for example ".in"
{

    char newfile[1000];
    snprintf(newfile,sizeof newfile,"%s_%s%s",pp,xx,zz);
    return(newfile);

}


void manager(int number_of_processors, int argc , char * argv[])
{
    int NMZ_MIN , NMZ_MAX , NB_MIN , NB_MAX , T_MIN , T_MAX, A_0_MIN,A_0_MAX;
    int NMZ_SIZE , NB_SIZE , T_SIZE, A_0_SIZE; //range of values of NMZ , NB, T
    int NMZ , T , NB, A_0, i, job[4]; //values to be passed to the eos function known as myfunction()

    i=number_of_processors;
    /* ASSIGNMENT */
    T_MIN=atoi(argv[1]);
    T_MAX=atoi(argv[2]);
    NB_MIN=atoi(argv[3]);
    NB_MAX=atoi(argv[4]);
    NMZ_MIN=atoi(argv[5]);
    NMZ_MAX=atoi(argv[6]);
    A_0_MIN=atoi(argv[7]);
    A_0_MAX=atoi(argv[8]);
    debug=atoi(argv[9]);

    /* determining ranges */
    NB_SIZE=(NB_MAX-NB_MIN)+1;
    NMZ_SIZE=(NMZ_MAX-NMZ_MIN)+1;
    T_SIZE=(T_MAX-T_MIN)+1;
    A_0_SIZE=(A_0_MAX-A_0_MIN)+1;

    int temp, who, tag;
    MPI_Status status;

    long int MAX_SIZE = (NMZ_SIZE)*(NB_SIZE)*(T_SIZE)*(A_0_SIZE);

    /*   for (i=0; i<number_of_processors; i++)
        {
            NMZ=(i%NMZ_SIZE)+NMZ_MIN;
            NB=((i/NMZ_SIZE)%NB_SIZE)+NB_MIN;
            T=((i/NB_SIZE)/NMZ_SIZE) + T_MIN;

            job[0] = NMZ;
            job[1]= NB;
            job[2] = T;

            MPI_Send(&job, 3, MPI_INT, i+1, i, MPI_COMM_WORLD);
        }
    */
    while (i<MAX_SIZE)
    {
        NMZ=(i%NMZ_SIZE)+NMZ_MIN;
        NB=((i/NMZ_SIZE)%NB_SIZE)+NB_MIN;
        T=((i/NB_SIZE)/NMZ_SIZE)%T_SIZE + T_MIN;
        A_0=(((i/NB_SIZE)/NMZ_SIZE)/T_SIZE) + A_0_MIN;
        if (!((NMZ<A_0) && (NMZ> (-1)*A_0)))
        {
            i++;
            continue;
        }
        else
        {
            job[0] = NMZ;
            job[1]= NB;
            job[2] = T;
            job[3] = A_0;

            MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            who = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            MPI_Send(&job, 4, MPI_INT, who, 1, MPI_COMM_WORLD);
            i++;
        }
    }

    for (i=0; i<number_of_processors; i++)
    {
        MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        who = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        MPI_Send(&job, 4, MPI_INT, who, 0, MPI_COMM_WORLD);
    }

    // MPI_Send(&temp,1, MPI_INT, 1, 1, MPI_COMM_WORLD);
}


void worker (int mpitask_id, long int MAX_SIZE, int task[])
{
    int result, tag,j;
    int NMZ , T , NB, A_0;
    MPI_Status status;
    int loopcount=0;
    int number_of_processors;

    double  initialtime, finaltime;
    double mytime;

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);

    int job[4];
    for(j=0; j<4; j++)
        job[j] = task[j];



    char buf1[10],buf2[10],buf3[10];
    char *temperature="T";
    char *NNBB="NB";
    char *NNMMZZ="NMZMAX";
    char *nnuucc="nuc";
    char *filetype1=".in";
    char dirname[100];
    char constfile[100];
    char nucfile[100];
    int ab=0;

    double *buf, *bbuf;
    int s1;

    int bufsize, bsize;

    /// MPI_Recv(&job, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
    A_0 = job[3];
    T = job[2];
    NB = job[1];
    NMZ = job[0];

    if (!((NMZ<A_0) && (NMZ > (-1)*A_0)))
    {
        MPI_Send(&result, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(&job, 4, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
        tag = status.MPI_TAG;
    }
    else tag = mpitask_id;
    initialtime=MPI_Wtime();


    int ipara,inuc,imb,inse,in_g,in_e,in_mu,in_tau,in_nu,
        in_cl[N_CL],in_hyp,in_tmes,in_ton,in_phase[N_PH],in_r,in_tab;

    XXX = 0.;

    /***** INPUT *****/

    /* RMF parametrization */
    /* 0: no interaction, 2: DD2, 3: TM1, 4: DD-MEdelta, default: DD */
    ipara = 2;

    /* calculation of nuclei or EoS */
    /* 1: nuclei */
    /* else: EoS*/
    inuc = 1;

    /* for calculation of EoS */
    /* 1: NSE calculation, default: RMF */
    /* noch gebraucht? in_ton reicht */
    inse = 0;

    /* for light clusters */
    /* 1: use non-relativistic Maxwell-Boltzmann statistics */
    /* 2: use relativistic Maxwell-Boltzmann statistics */
    /* else: use original statistics */
    imb = 1;

    /* include photon */
    in_g = 0;
    /* include electron */
    in_e = 0;
    /* include muon */
    in_mu = 0;
    /* include tauon */
    in_tau = 0;
    /* include neutrinos */
    in_nu = 0;
    /* include deuteron */
    in_cl[0] = 0;
    /* include triton */
    in_cl[1] = 0;
    /* include helion */
    in_cl[2] = 0;
    /* include alpha */
    in_cl[3] = 0;
    /* include np 3S1 channel */
    in_cl[4] = 0;
    /* include np 1S0 channel */
    in_cl[5] = 0;
    /* include nn 1S0 channel */
    in_cl[6] = 0;
    /* include pp 1S0 channel */
    in_cl[7] = 0;
    /* include triton continuum */
    in_cl[8] = 0;
    /* include helion continuum */
    in_cl[9] = 0;
    /* include alpha continuum */
    in_cl[10] = 0;
    /* include hyperons */
    in_hyp = 0;
    /* include thermal mesons */
    in_tmes = 0;
    /* include table of nuclei (only for eos calculation) */
    in_ton = 0;
    if (inuc == 1) in_ton = 0;

    /* include drop phase */
    in_phase[1] = 0;
    /* include bubble phase */
    in_phase[2] = 0;
    /* include hole phase */
    in_phase[3] = 0;

    /* 1: read data2files for radii */
    in_r = 0;

    /* 1: standard input file eos_table.in
       else: input file eos_tab.in */
    in_tab = 1;

    /* 0: dependence of energy shifts on meson fields
       else: dependence of energy shifts on densities */
    DBDEP = 0;
    int kk=0;
    /***** INITIALIZATION *****/

    snprintf(constfile,sizeof constfile,"%s%s_%s%s_%s%s_%d",temperature,buf1,NNBB,buf2,NNMMZZ,buf3,mpitask_id);

    init_eos(ipara,imb,in_ton,constfile);

    int old_n_b = job[1];
    double nb_offset, offset;
    nb_offset = (mpitask_id-1)*sizeof(char)*150;
    offset = (mpitask_id-1)*sizeof(char)*150;

    while (tag)
    {
        //mytime=MPI_Wtime();

        A_0 = job[3];
        T = job[2];
        NB = job[1];
        NMZ = job[0];

        snprintf(buf1,10,"%d",T);
        snprintf(buf2,10,"%d",NB);
        snprintf(buf3,10,"%d",NMZ);

        snprintf(constfile,sizeof constfile,"%s%s_%s%s_%s%s_%d",temperature,buf1,NNBB,buf2,NNMMZZ,buf3,mpitask_id);
        // snprintf(nucfile,sizeof(nucfile),"T_%d_nuc.dat",T);
        //job.out files created for each task
        //myfile=fopen(changestring("J",constfile,".out"),"a+");

        //CALL CODE
        data2 * point;
        point=(data2*)malloc(sizeof(data2));
        //point=myfunction(constfile,NMZ,NB,T,A_0,point);


        /***** CALCULATION *****/

        if (inuc == 1)
        {
            point=get_nuc_table(inse,in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,NMZ,NB,T,constfile,point, A_0);
        }
        else
        {
            get_eos_table(inse,in_g,in_e,in_mu,in_tau,in_nu,in_cl,in_hyp,in_tmes,
                          in_ton,in_phase,in_r,in_tab,constfile);
        }

        //fclose(myfile);


        //    MPIO_Request *req = (MPIO_Request*)malloc(count * sizeof(MPIO_Request));
        //   MPIO_Request *req1 = (MPIO_Request*)malloc(count * sizeof(MPIO_Request));


        /*    mytime=MPI_Wtime()-mytime;

            FILE *sakfile;
            sakfile=fopen("mytime.txt","a+");
            fprintf(sakfile,"%d %d %d %d %lf %d\n",NMZ,T,NB,A_0,mytime,mpitask_id);
            fclose(sakfile);
          */
        int reqcount = count;

        MPIO_Request *req = (MPIO_Request*)malloc(reqcount * sizeof(MPIO_Request));
        MPIO_Request *req1 = (MPIO_Request*)malloc(reqcount * sizeof(MPIO_Request));

        MPI_File fh, fh1;

        MPI_Info inf;

        MPI_Info_create(&inf);
        MPI_Info_set(inf, "noncoll_write_bufsize", "300000");

        sprintf(nucfile,"T_%i.dat", T);



        for(loopcount=0; loopcount<count; loopcount++)
        {
            char tstring[150];
            char atfile[50];

            int write_count = sprintf(tstring, "\n%e %e %i %i %e %e %e %e %e %e %e %e\n", point->a[loopcount],point->b[loopcount],point->c[loopcount],point->d[loopcount],point->e[loopcount],point->f[loopcount],point->g[loopcount],point->h[loopcount],point->i[loopcount],point->j[loopcount],point->k[loopcount],point->l[loopcount]);

            if(old_n_b != NB)
            {
                nb_offset = (mpitask_id-1)*sizeof(char)*150;
                old_n_b=NB;
            }

            sprintf(atfile,"NB_%i.dat", NB);

            MPI_File_open(MPI_COMM_SELF,nucfile,(MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND ), inf, &fh);


            MPI_File_iwrite_at( fh, offset ,tstring, write_count, MPI_CHAR, &req[loopcount]);

            MPI_File_close( &fh );
            offset += number_of_processors*sizeof(char)*150;
            MPI_File_open(MPI_COMM_SELF,atfile,(MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND), inf, &fh1);

            MPI_File_iwrite_at(  fh1, offset ,tstring, write_count, MPI_CHAR,&req1[loopcount]);

            MPI_File_close( &fh1 );

            nb_offset += number_of_processors*sizeof(char)*150;
        }

        free(point);




        // MPI_Waitall(count, req, MPI_STATUSES_IGNORE);
        //      MPI_Waitall(count, req1, MPI_STATUSES_IGNORE);
        resetcount();
        //free(req);




        MPI_Send(&result, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Recv(&job, 4, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
        tag = status.MPI_TAG;

        MPI_Waitall(reqcount, req, MPI_STATUSES_IGNORE);
              MPI_Waitall(reqcount, req1, MPI_STATUSES_IGNORE);

        free(req);
        free(req1);
//MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
    /*    finaltime=MPI_Wtime()-initialtime;
        FILE *a;

        a=fopen("totaltime.txt","a+");
        fprintf(a,"%d %lf\n",mpitask_id,finaltime);
        fclose(a);
       MPI_Buffer_detach(&bbuf,&bsize);
    */

}
