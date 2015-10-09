#include "Variable.h"


class Get_k2_red
{
public:
	Get_k2_red();
	~Get_k2_red();
	double get_k0_red(double x);
	double get_k1_red(double x);
	double get_k2_red(double x);
};

class Worker
{
/*private:
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
	double XXX;
	int DBDEP;
	*/
public:
	Worker();
	std::vector <double> worker(int mpitask_id, int NMZ_MIN, int NMZ_MAX, int NB_MIN, int NB_MAX, int T_MIN, int T_MAX, int T_DISP, int A_0_DISP, int NB_DISP, int NMZ_DISP, int A_0_MIN, int A_0_MAX);
	data2* myfunction(char *s, int t_nmz, int t_nb, int t_t, int a_0, data2*point);
	data2* work(int task[], int size,int mpitask_id, int displacement,double *worked_for,int displacement1);
	~Worker();

	friend class Slave;
};


class Fold_dens
{
public:
	Fold_dens();
	~Fold_dens();
	int fold_dens(int ip);
	double func_fold(double r, double rp, double a2);
};

class Solve_rmf
{

public:
	Solve_rmf();
	~Solve_rmf();
	int solve_rmf(int iph, double r);
	int init_discr_ws_cell(int iph);
	int init_dens(void);
	double f_ws(double r, double rr, double a);
	int init_mf(void);
	int get_self_energies(int ic);
	int init_self_energies(void);
	int get_anz_nuc(void);

	friend int get_nucleus(double a, double z, char *cnsfile);
	friend 	int get_eos_composition(int in_cell, char* cnsfile);
	friend int fit_radius2(int iph, double r_in);

};

class Slave
{
public:
	Slave();
	void slave(data3 *table, int NMZ_MIN, int NMZ_MAX, int NB_MIN, int NB_MAX, int T_MIN, int T_MAX, int T_DISP, int A_0_DISP, int NB_DISP, int NMZ_DISP, int A_0_MIN, int A_0_MAX);
	~Slave();
};

class Master
{

public:
	Master();
	void master(data3 * table, int a_size);
	~Master();
};

class Manager
{
public:
	Manager();
	friend class Slave;
	void manager(int number_of_processors, int NMZ_MIN, int NMZ_MAX, int NB_MIN, int NB_MAX, int T_MIN, int T_MAX, int T_DISP, int A_0_DISP, int NB_DISP, int NMZ_DISP, int A_0_MIN, int A_0_MAX);
	~Manager();
};

class Init_eos
{

public:
	Init_eos();
	~Init_eos();
	int init_parameters(void);
	int init_eos(int ipara, int imb, int in_ton, char*cnsfile);
	int init_rmf(int ipara);
	int init_particles(int imb);
	int init_ton(int in_ton, char *cnsfile);
	int init_nuclei_ame03(int in_ton, char *cnsfile);
	friend data2* myfunction(char *s, int t_nmz, int t_nb, int t_t, int a_0, data2*point);
};

class Init_be_cl
{
public:
	Init_be_cl();
	~Init_be_cl();
	int init_be_cl(void);
	double experfc(double z);
};

class Get_rho_part :
	public Get_k2_red
{
public:
	Get_rho_part();
	~Get_rho_part();
	double get_rho_part(double mu, int ip, int ic);
	int get_rho(double mu);
	int get_rhos(double mu);
	double get_simint(int order_max, double prec, int ifunc,
		double x_min, double x_max, double rpar, int ipar);
	int get_simint_bounds_f(double mu, int ipa);
	int get_simint_bounds_b(double mu, int ipa);
	double inv_disp_rel(double e, int ipa);
	double disp_rel(double k, int ipa);

	friend double get_n_lep(double mu, double fc, int ic);
	friend int get_nb_yq(double mu_b, double mu_q, int ic);
};

class Get_properties :
	public Get_k2_red
{
public:
	Get_properties();
	~Get_properties();
	int get_properties(void);
	int get_prop_part(int ip);
	int get_prop_ton(void);
	int get_prop(double mu);

	friend int solve_rmf(int iph, double r);
};

class Get_nucleus
{

public:
	Get_nucleus();
	~Get_nucleus();
	int get_nucleus(double a, double z, char *cnsfile);
	int init_sol(void);
	int iexs_sol(int iph);
	int save_sol(int iph);

	friend int get_eos_local(double t, double n_b, double y_q, int inse, int in_g, int in_e, int in_mu, int in_tau, int in_nu, int *in_cl, int in_hyp, int in_tmes, int in_cell, int in_ton, int inuc, double a, double z, char *cnsfile);
	friend int get_eos_composition(int in_cell, char* cnsfile);
};

class Get_nuc_table : public Init_be_cl
{


public:
	Get_nuc_table();
	~Get_nuc_table();
	data2* get_nuc_table(int inse, int in_g, int in_e, int in_mu, int in_tau, int in_nu,
		int *in_cl, int in_hyp, int in_tmes, int t_nmz, int t_nb, int t_t, char *cnsfile, data2* point, int a_0);
	int get_eos_local(double t, double n_b, double y_q, int inse, int in_g, int in_e, int in_mu, int in_tau, int in_nu, int *in_cl, int in_hyp, int in_tmes, int in_cell, int in_ton, int inuc, double a, double z, char *cnsfile);
	int get_results_ws_cell(char *cnsfile);

	friend data2* myfunction(char *s, int t_nmz, int t_nb, int t_t, int a_0, data2*point);
	friend int get_nucleus(double a, double z, char *cnsfile);
};

class Get_mf :
	public Fold_dens
{
public:
	Get_mf();
	~Get_mf();
	int get_mf(int ic);
	int get_mf_greens(int imf, double *tmps);
	int get_acc_mf(int imf);
	double simpson_mod(int np, double *vec);
	int get_sources(void);
	int init_acc_mf(int imf);
	int get_ton_g(int ip);
	int get_ton_mass(double n_p, double n_n);
	int get_ton_mass_shift(int ip, double n_ref, double n_eff);

	friend int get_self_energies(int ic);
};

class Get_mat_inv
{
public:
	Get_mat_inv();
	~Get_mat_inv();
	int get_mat_inv(int dim);
	int improve(int n);
	int lu_backsubstitution(int n);
	int lu_decomposition(int n);

	friend 	int get_acc_mf(int imf);
};

class Get_eos_table : public Get_nuc_table
{

public:
	Get_eos_table();
	~Get_eos_table();
	int get_eos_table(int inse, int in_g, int in_e, int in_mu, int in_tau, int in_nu,
		int *in_cl, int in_hyp, int in_tmes, int in_ton,
		int *in_phase, int in_r, int in_tab, char *cnsfile);
	int get_parameter_tables(int ir_y, int ir_n, int ir_t,
		double d_y_q, double f_n_b, double t_ref, double f_t, char *cnsfile);
	double get_eos_local_beta(double t, double n_b, int inse,
		int in_g, int in_e, int in_mu, int in_tau, int in_nu,
		int *in_cl, int in_hyp, int in_tmes,
		int in_cell, int in_ton, int inuc, double a, double z, char* cnsfile);
	friend data2* myfunction(char *s, int t_nmz, int t_nb, int t_t, int a_0, data2*point);
};

class Get_eos_composition :
	public Fold_dens
{

public:
	Get_eos_composition();
	~Get_eos_composition();
	int get_eos_composition(int in_cell, char* cnsfile);
	int get_fractions(void);
	int get_coul_corr(void);
	int select_sol(int ip_max);
	friend int get_eos_local(double t, double n_b, double y_q, int inse, int in_g, int in_e, int in_mu, int in_tau, int in_nu, int *in_cl, int in_hyp, int in_tmes, int in_cell, int in_ton, int inuc, double a, double z, char *cnsfile);
};

class Get_chemical_potential_lep
{
public:
	Get_chemical_potential_lep();
	~Get_chemical_potential_lep();
	double get_n_lep(double mu, double fc, int ic);
	double get_chemical_potential_lep(void);

	friend int solve_rmf(int iph, double r);
};

class Get_chemical_potential_bar
{
public:
	Get_chemical_potential_bar();
	~Get_chemical_potential_bar();
	int get_chemical_potential_bar(double mu_b, double mu_q);
	int fit_asymmetry(double mu_b, double dmu_b, double mu_q, double dmu_q);
	int fit_density(double mu_b, double dmu_b, double mu_q);
	int get_nb_yq(double mu_b, double mu_q, int ic);
	int get_dens_ton(int ic);
	int heapsort(int n);

	friend int solve_rmf(int iph, double r);
};

class Get_be_cl3
{
public:
	Get_be_cl3();
	~Get_be_cl3();
	int get_be_cl3(int i1, double n_p, double n_n);
	int get_dbe_p_g(int icl);
	int get_dbe_p_j(void);
	int get_shift(double de, double det, double x, double x0);

	friend int get_self_energies(int ic);
	friend int get_properties(void);
};

class Fit_radius2
{
public:
	Fit_radius2();
	~Fit_radius2();
	int fit_radius2(int iph, double r_in);
	double solve_a(int ic, int iph, double a);
	int fit_aa4_new(int iph, double r);

	friend 	int get_eos_composition(int in_cell, char* cnsfile);
};

class F_func
{
public:
	F_func();
	~F_func();
	double f_func(double k, int ifunc, double mu, int ipa);
	double func_fermi(double k, double mu, int ipa);
	double func_bose(double k, double mu, int ipa);

	friend 	double get_simint(int order_max, double prec, int ifunc,
		double x_min, double x_max, double rpar, int ipar);
};

class Cycle_solve :
	public Fold_dens, public Get_mat_inv
{
public:
	Cycle_solve();
	~Cycle_solve();
	int cycle_solve(int iph, double mu_b, double mu_q);
	int get_acc(void);
	int get_acc0(void);
	int init_acc(void);
	int init_acc0(void);

	friend int solve_rmf(int iph, double r);
};

class Main_class
{
public:
	Main_class();
	~Main_class();
	int main(int argc, char * argv[]);
};

