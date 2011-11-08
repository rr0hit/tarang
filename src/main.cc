/*! \file  main.cc
 *
 *	@brief  Main file for Tarang (sequential). Calls various modules.
 *  @sa main.h
 *
 * @note Reads which turbulence to start (Fluid, MHD, Passive scalcar, RBslip),
 *					the dimensionality of simulation (2D or 3D), and
 *					which directory to use for input and output.
 *
 *	@author  M. K. Verma
 *	@version 4.0 MPI
 *	@date Sept 2008
 */


#include "main.h"

// declarations of global MPI,  fftw_plan and random variables

int my_id;								// My process id
int numprocs;							// No of processors
const int master_id = 0;				// Id of master proc
ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
ptrdiff_t local_N2, local_N2_start;
MPI_Status status;

// for fftw_original
fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;

// for split fftw
fftw_plan r2c_2d_plan_FOUR, c2r_2d_plan_FOUR;
fftw_plan c2c_1d_forward_plan_FOUR, c2c_1d_inverse_plan_FOUR;

// for SCFT
fftw_plan r2c_plan_SCFT, c2r_plan_SCFT, sintr_plan_SCFT, costr_plan_SCFT,
	isintr_plan_SCFT, icostr_plan_SCFT;  // i for inverse

fftw_plan r2c_1d_plan_SCFT, c2r_1d_plan_SCFT;  // for 2D

int		globalvar_fftw_original_switch = 0;
int		globalvar_anisotropy_switch;			// 1,2,3 for x,y,z directions
int		globalvar_waveno_switch;				// 0 for actual (default), 1 for grid


string globalvar_prog_kind;

// for RB
string	globalvar_Pr_switch;					// Prandtl number switch (PRLARGNE..) for RB
string	globalvar_RB_Uscaling;					// UBscaling (ULARGE... ) for RB
DP		globalvar_Ra;							// Rayleigh number
DP		globalvar_r;							// normalized Rayleigh number
DP		globalvar_Pr;							// Prandtl number
DP		globalvar_temperature_grad;				// +1 for convection; -1 for stratification;
// factor for u3 in temperature eqn
DP		globalvar_alpha_DT;						// nondim number alpha*DT

Uniform<DP> SPECrand;					// Global variable for random no generation


int main(int argc, char** argv)
{

	string		data_dir_name;
	time_t		start;
	time_t		end;

	SPECrand.seed((unsigned int)time(0));				// Initialize random number seed

	// MPI INIT
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	time (&start);
	Read_prog_para(globalvar_prog_kind, data_dir_name);

	if (globalvar_prog_kind == "INC_FLUID")
		Ifluid_main(data_dir_name);
	else if (globalvar_prog_kind == "INC_FLUID_DIAG")
		Ifluid_diag_main(data_dir_name);
	else if (globalvar_prog_kind == "INC_SCALAR")
		Iscalar_main(data_dir_name);
	else if (globalvar_prog_kind == "INC_SCALAR_DIAG")
		Iscalar_diag_main(data_dir_name);
	else if (globalvar_prog_kind == "INC_MHD")
		IMHD_main(data_dir_name);
	else if (globalvar_prog_kind == "INC_MHD_DIAG")
		IMHD_diag_main(data_dir_name);
	else if (globalvar_prog_kind == "RB_SLIP")
		RB_slip_main(data_dir_name);
	else if (globalvar_prog_kind == "RB_SLIP_DIAG")
		RB_slip_diag_main(data_dir_name);
	else if (globalvar_prog_kind == "RB_SLIP_MHD")
		RB_slipMHD_main(data_dir_name);
	else if (globalvar_prog_kind == "NON_BOUSSINESQ")
		NonBoussinesq_main(data_dir_name);
	else
		if (my_id == master_id)
			cout << "ERROR: NO CHOICE (e.g., INC_FLUID) WAS MADE " << endl;

	time(&end);
	double dif = difftime(end, start);
	if (my_id == master_id)
		cout << endl << "PROGRAM TERMINATING HERE" << endl << endl
		     << "TOTAL TIME ELAPSED: " << dif << " sec" << endl;

	MPI_Finalize();
}
