#include "../main.h"

// Global variables
extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;		// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2, local_N2_start;
extern MPI_Status status;

extern Uniform<DP> SPECrand;					// for random vars--  declared in main.cc


//*********************************************************************************************

int RB_slip_diag_main(string data_dir_name)
{

	cout << "ENTERING RB_slip_main :   my_id " << my_id << endl;

	//**************** Variable declarations *******************

	// Variables connected to the field variables and integration

	int			N[4];
	int			no_of_fields = 2;				// U, B fields
	DP			diss_coefficients[2];
	DP			hyper_diss_coefficients[2];


	ifstream	para_file;					// field and input-output parameters stored here

	string		string_switches[MAXSIZE_STRING_SWITCHES_ARRAY];

	Array<int,1>	switches( MAXSIZE_SWITCHES_ARRAY);

	// solver para
	Array<int,1>	solver_meta_para(MAXSIZE_SOLVER_META_PARA);
	Array<int,1>	solver_int_para(MAXSIZE_SOLVER_INT_PARA);
	Array<DP,1>		solver_double_para(MAXSIZE_SOLVER_DOUBLE_PARA);
	string			solver_string_para[MAXSIZE_SOLVER_STRING_PARA];


	Array<int,1>	diagnostic_procedure(MAXSIZE_DIAGNOSTIC_PARA);
	diagnostic_procedure = 0;

	// Variables connected to the output functions

	//! Tinit, Tfinal, Tdt, Tdiagnostic_int in index 1..4
	Array<DP, 1>time_para(MAXSIZE_TIME_PARA_ARRAY);
	time_para = 0.0;

	Array<DP,1>		time_save(MAXSIZE_TIME_SAVE_ARRAY);
	time_save		= 0.0;

	Array<int,1>	misc_output_para(MAXSIZE_MISC_OUTPUT_PARA);

	//!  contains no_output_k, no_output_r
	Array<int,1>	no_output_k_r(3);
	Array<int, 2>	out_k_r_array(MAXSIZE_OUT_K_R_ARRAY,4);
	out_k_r_array		= 0;			// intialize



	// Variables connected to the energy transfer

	Array<int,1>	ET_parameters(MAXSIZE_ET_PARAMETERS);

	Array<DP,1>		ET_shell_radii_sector_array(MAXSIZE_ET_RADII_SECTOR_ARRAY);

	ET_shell_radii_sector_array = 0.0;



	// init-cond para
	Array<int,1>	init_cond_meta_para(MAXSIZE_INIT_COND_META_PARA);
	Array<int,1>	init_cond_int_para(MAXSIZE_INIT_COND_INT_PARA);
	Array<DP,1>		init_cond_double_para(MAXSIZE_INIT_COND_DOUBLE_PARA);
	string			init_cond_string_para[MAXSIZE_INIT_COND_STRING_PARA];


	// forcing para
	Array<int,1>	force_meta_para(MAXSIZE_FORCE_META_PARA);
	Array<int,1>	force_int_para(MAXSIZE_FORCE_INT_PARA);
	Array<DP,1>		force_double_para(MAXSIZE_FORCE_DOUBLE_PARA);
	string			force_string_para[MAXSIZE_FORCE_STRING_PARA];


	//**************** Program starts here ****************

	// Read field para and input-output parameters

	string		filename = "/in/para.d";
	filename = data_dir_name+ filename;
	para_file.open(filename.c_str());		// filename = data_dir_name/in/para.d

	Read_para(para_file, 3, no_of_fields, N, string_switches, switches,
		  diss_coefficients, hyper_diss_coefficients,
		  solver_meta_para, solver_int_para, solver_double_para, solver_string_para,
		  time_para, time_save,
		  misc_output_para, ET_parameters, ET_shell_radii_sector_array,
		  no_output_k_r, out_k_r_array,
		  init_cond_meta_para, init_cond_int_para,
		  init_cond_double_para, init_cond_string_para,
		  force_meta_para, force_int_para, force_double_para, force_string_para);

	string_switches[0] = data_dir_name;

	globalvar_anisotropy_switch = switches(14);
	globalvar_waveno_switch = switches(15);



	if (my_id == master_id)
		cout << " ================ RB PARAMTERS ==================== "  << endl;

	// kfactor computation
	DP	k0 = M_PI/sqrt(2.0);

	DP	kfactor[4];

	if (solver_int_para(1) == 0) {	// kfactor as given
		kfactor[1] = solver_double_para(1);
		kfactor[2] = solver_double_para(2);
		kfactor[3] = solver_double_para(3);

	} else if (solver_int_para(1) == 1) {	// kfactor multiplied by some predefined constant
		kfactor[1] = M_PI * solver_double_para(1);
		kfactor[2] = k0 * solver_double_para(2);
		kfactor[3] = k0 * solver_double_para(3);
	}

	if (my_id == master_id)
		cout << "kfactor: " << kfactor[1] << " "  << kfactor[2] << " "  << kfactor[3]
		     << endl << endl;

	// Ra, Pr, Pr_switch,
	if (solver_int_para(2) == 0)
		globalvar_Ra = solver_double_para(4);
	else if (solver_int_para(2) == 1)
		globalvar_Ra = (27.0*pow4(M_PI)/4.0)*solver_double_para(4);

	globalvar_Pr = solver_double_para(5);

	globalvar_temperature_grad = solver_double_para(6);  // +1 for RB and -1 for stratified

	if (my_id == master_id) {
		cout << " Rayleigh no:  Ra =  " << globalvar_Ra << endl;
		cout << " Prandtl no:  Pr = " << globalvar_Pr << endl;
		cout << " Temperature gradient (+1 for RB, -1 for stratififed flow) = "
		     << globalvar_temperature_grad << endl;
	}

	// Pr_switch and RB_Uscaling
	globalvar_Pr_switch = solver_string_para[1];

	globalvar_RB_Uscaling = solver_string_para[2];

	if (my_id == master_id) {
		cout << " Pr_switch = " << globalvar_Pr_switch << endl;
		cout << " RB_Uscaling = " << globalvar_RB_Uscaling << endl;
		cout << "====================================================" << endl << endl;
	}

	// Construct output prefix for all the output files
	string prefix_str,  Pr_str, r_str;
	ostringstream Pr_buffer, r_buffer;

	Pr_buffer << globalvar_Pr;
	Pr_str = Pr_buffer.str();

	r_buffer << globalvar_Ra;
	r_str = r_buffer.str();

	prefix_str = "%% Pr = " + Pr_str +  "  r = " + r_str;	// goes into all the output files


	// Set Dissipation coefficients
	if (globalvar_Pr_switch == "PRZERO") {
		diss_coefficients[0] = 1.0;
		diss_coefficients[1] = 0.0;
	}

	else if (globalvar_Pr_switch == "PRLARGE") {
		if (globalvar_RB_Uscaling == "USMALL") {
			diss_coefficients[0] = globalvar_Pr;              //  Coeff of grad^2 u
			diss_coefficients[1]  = 1.0;			// Coeff of grad^2 T
		}
		else if (globalvar_RB_Uscaling == "ULARGE") {
			diss_coefficients[0] = sqrt(globalvar_Pr/globalvar_Ra);
			diss_coefficients[1]  = 1/sqrt(globalvar_Pr*globalvar_Ra);
		}
	}

	else if (globalvar_Pr_switch == "PRSMALL") {
		if (globalvar_RB_Uscaling == "USMALL") {
			diss_coefficients[0] = 1.0;
			diss_coefficients[1]  = 1/globalvar_Pr;
		} else if (globalvar_RB_Uscaling == "ULARGE") {
			diss_coefficients[0] = sqrt(globalvar_Pr/globalvar_Ra);
			diss_coefficients[1]  = 1/sqrt(globalvar_Pr*globalvar_Ra);
		}
	}

	else if (globalvar_Pr_switch == "PRINFTY") {
		if (globalvar_RB_Uscaling == "USMALL") {
			diss_coefficients[0] = globalvar_Pr;              //  Coeff of grad^2 u
			diss_coefficients[1]  = 1.0;			// Coeff of grad^2 T
		} else if (globalvar_RB_Uscaling == "ULARGE") {
			cout << "ERROR: For PRINFTY, ULARGE option is not allowed" << endl;
			exit(0);
		}
	}





	string basis_type  = string_switches[1];

	// Local_N1, local_N2 assignments
	if (N[2] > 1) {
		if ((basis_type == "FOUR") && (globalvar_fftw_original_switch == 1))	{
			int alloc_local;
			alloc_local = fftw_mpi_local_size_3d_transposed(N[1], N[2], N[3], MPI_COMM_WORLD,
									&local_N1, &local_N1_start, &local_N2, &local_N2_start);
		} else {
			local_N1 = N[1]/numprocs;
			local_N2 = N[2]/numprocs;
			local_N1_start = my_id * local_N1;
			local_N2_start = my_id * local_N2;
		}
	}

	else if (N[2] == 1) {   // Work on it
		if (basis_type == "FOUR") {
			int alloc_local;
			alloc_local = fftw_mpi_local_size_2d_transposed(N[1], N[3], MPI_COMM_WORLD,
									&local_N1, &local_N1_start, &local_N2, &local_N2_start);
		} else if (basis_type == "SCFT") {
			local_N1 = N[1]/numprocs;
			local_N2 = N[3]/numprocs;
			local_N1_start = my_id * local_N1;
			local_N2_start = my_id * local_N2;
		}
	}

	if (my_id == master_id) cout << "No or processors = " << numprocs << endl << endl;

	cout << "MY ID, local_N1, local_N1_start, local_N2, local_N2_start = " << my_id << "  "
	     << local_N1 << "  " << local_N1_start << "  " << local_N2 << "  " << local_N2_start
	     << endl << endl;



	// Constructors of Vector and Scalar fields

	IncFluid  U(N, string_switches, switches, kfactor,
		    diss_coefficients[0], hyper_diss_coefficients[0],
		    solver_meta_para, solver_int_para, solver_double_para, solver_string_para,
		    time_para, time_save,
		    misc_output_para, no_output_k_r, out_k_r_array,
		    ET_parameters, ET_shell_radii_sector_array,
		    init_cond_meta_para, init_cond_int_para,
		    init_cond_double_para, init_cond_string_para,
		    force_meta_para, force_int_para, force_double_para, force_string_para
		    );

	IncSF     T(N, string_switches, switches, kfactor,
		    diss_coefficients[1], hyper_diss_coefficients[1],
		    misc_output_para
		    );


	// Create FFTW plans

	U.Init_fftw_plan();

	// Initialize vector fields;

	if (my_id == master_id)
		U.Open_field_input_files();

	U.Read_init_cond(T);

	if (my_id == master_id)
		U.Close_field_input_files();

	//
	//

	if (my_id == master_id)	{
		U.Open_output_files(prefix_str);
		cout << endl << "STARTING THE SIMULATION NOW" << endl;
	}


	U.Tnow = U.Tinit;


	int i=1;
	while (diagnostic_procedure(i) < MY_MAX_INT) {
		switch (diagnostic_procedure(i)) {
		case (0) : U.Output_global(T);		break;
		case (1) : U.Output_shell_spectrum(T); break;
		case (2) : U.Output_ring_spectrum(T);  break;
		case (3) : U.Output_cylinder_ring_spectrum(T);  break;
		case (4) : U.Output_flux(T);		break;
		case (5) : U.Output_shell_to_shell(T); break;
		case (6) : U.Output_cylinder_ring_to_ring(T);  break;
		case (7) : U.Output_structure_fn(T);  break;
		case (8) : U.Output_planar_structure_fn(T);  break;
		case (9) : U.Output_moment(T);  break;
		case (10) : U.Output_field_k(T);  break;
		case (11) : U.Output_field_r(T);  break;
		case (12) : U.Output_Skpq(T);  break;
		case (13) : U.Output_realfield(T);  break;
		case (14) : U.Output_field_reduced(T); break;
		}

		i++;
	}


	if (my_id == master_id)
		U.Close_output_files();

	return 1;
}
