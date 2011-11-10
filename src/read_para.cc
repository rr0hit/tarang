#include "main.h"

extern int my_id;								// My process id
extern int numprocs;							// No of processors
extern const int master_id;						// Id of master proc
extern ptrdiff_t local_N1, local_N1_start;			// N1 size and start of i1 in the currentproc
extern ptrdiff_t local_N2;
extern MPI_Status status;

//*********************************************************************************************

/** @brief Reads parameters
 *
 * @param  para_file	The name of the paramter file (dir/in/para.d)
 * @param  dim	Dimension of the field array (2, 3)
 * @param  number_of_fields	Number of fields in the simulation (e.g., fluid: 1, MHD: 2)
 * @param  N[] The size of field arrays (real space)
 *
 * @param string_switches[1] Basis_type_switch (FOUR, SCFT)
 * @param string_switches[2] alias_switch (ALIAS, DEALIAS)
 * @param string_switches[3] Integration scheme (EULER, RK2, RK4)
 * @param string_switches[4] Input_number_mode (ASCII, BINARY)
 * @param string_switches[5] Output_number_mode (ASCII, BINARY).
 *
 * @param switches(1) Hyper_dissipation switch (0: no, 1: yes)
 * @param switches(2) Low_dimensional switch
 * @param switches(3) LES switch
 * @param switches(4) Anisotropy_ring switch
 * @param switches(5) Anisotropy_cylinder switch
 * @param swithces(6) structure_fn_switch
 * @param switches(7) planar_structure_fn_switch
 * @param switches(8) skpq_switch
 * @param switches(9) output_real_space_switch:  Output V(r) at select points
 * @param switches(10) 1 if the horizontal basis function is cos(n k0 y) or sin(n k0 y).
 *						We choose Vx = sin(m pi x) cos(n k0 y)
 *
 * @param diss_coefficient[] dissipation coefficient for each field of number_of_fields.
 * @param hyper_diss_coefficient[]  Hyper dissipation coefficient for each field of
 *						number_of_fields. (keep it zero if the hyper_dissipation switch is 0).
 *
 * @param time_para(1) Initial time
 * @param time_para(2) Final time
 * @param time_para(3) dt_fixed (dt is minimum of dt_CFL & dt_fixed)
 * @param time_para(4) Time from which diagnostic functions will be evoked.
 * @param time_para(5) TIme from which Structure function calculation will be evoked.
 *
 *
 * @param time_save(1)	global variable writing interval (eg. total energy)
 * @param time_save(2)  field variable writing interval.
 * @param time_save(3)  Interval for writing field frequent (Overwritten- contains the most
 *							frequent field. Useful for powercuts and crashes)
 * @param time_save(4)  Interval for writing field variable in reduced grid.
 * @param time_save(5)  real field variable writing interval.
 * @param time_save(6)  Interval for writing field variables at select modes.
 * @param time_save(7)  Interval for writing real field variables at specified space locations.
 * @param time_save(8)  Interval for writing energy spectrum
 * @param time_save(9)  Interval for writing pressure spectrum
 * @param time_save(10) Interval for writing isotropic flux
 * @param time_save(11) Interval for writing shell-to-shell energy transfer
 * @param time_save(12) Interval for writing ring-to-ring spectrum.
 * @param time_save(13) Interval for writing ring-to-ring transfers.
 * @param time_save(14) Interval for writing cylinderical ring-to-ring spectrum.
 * @param time_save(15) Interval for writing cylinderical ring-to-ring transfers.
 * @param time_save(16)	Structure function calculation time interval.
 * @param time_save(17) Interval for writing moments of fields.
 * @param time_save(18) Interval for writing Skpq.
 * @param time_save(19) Interval for writing some global variables in cout.
 *
 * @param misc_output_para(1) Ring spectrum: sector input scheme (0: uniform angle,
 *											1: angle s.t. each sector has equal modes).
 * @param misc_output_para(2) Number of sectors for ring spectrum.
 * @param misc_output_para(3) Number of slabs in cylinders.
 * @param misc_output_para(4) Structure function q-min.
 * @param misc_output_para(5) Structure function q-max.
 * @param misc_output_para(6:5+dim) N_out_reduced[]
 *
 * @param ET_parameters(1) Energy transfer: real-imag-switch (0: no, 1: yes).
 * @param ET_parameters(2) Isotropic energy transfer: no of spheres.
 * @param ET_parameters(3) Isotropic energy transfer: no of shells.
 * @param ET_parameters(4) Ring energy transfer: no of ring-shells.
 * @param ET_parameters(5) Ring energy transfer: no of sectors.
 * @param ET_parameters(6) Cylindrical ring energy transfer: no of cylinder-shells.
 * @param ET_parameters(7) Cylindrical ring energy transfer: no of cylinder-slabs.
 * @param ET_parameters(8) ET Input scheme: shells
 * @param ET_parameters(9) ET Input scheme: ring-shells
 * @param ET_parameters(10) ET Input scheme: sector
 * @param ET_parameters(11) ET Input scheme: cylinder-shells.
 * @param ET_parameters(12) ET Input scheme: cylinder-kpll.
 *
 * @param ET_shell_radii_sector_array() shell radii, ring-shell-radii, cylinder-shell-radii,
 *											sectors if input scheme is so.
 *
 * @param no_output_k_r(1) Number of fourier-space points for which field values
 *												and Tk are to be printed.
 * @param no_output_k_r(2)  No of real-space points for which the field are to be printed.
 *
 * @param out_k_r_array()  range (1:no_output_k_r(1), 1:dim)  k values (grid wavenumber).
 * @param out_k_r_array()  range (no_output_k_r(1)+1:no_output_k_r(1)+no_output_k_r(2), 1:dim)
 *											r values.
 *
 * @param field_input_meta_para(1) Field input procedure.
 * @param field_input_meta_para(2) Number of init cond parameters
 * @param field_input_meta_para(3:2+dim) N_in_reduced[].
 *
 * @param init_cond_para() range(1:field_input_meta_para(2)) init cond parameters
 *
 * @param force_field_meta_para(1) Force field procedure
 * @param force_field_meta_para(2) Number of force field parameters.
 *
 * @param force_field_para() range(1:force_field_meta_para(2)) Force field parameters.
 *
 */

void Read_para
(
 ifstream& para_file,
 int dim,
 int number_of_fields,
 int N[],
 string string_switches[],
 Array<int,1> switches,
 double diss_coefficient[],
 double hyper_diss_coefficient[],
 Array<int,1> solver_meta_para,
 Array<int,1> solver_int_para,
 Array<DP,1> solver_double_para,
 string solver_string_para[],
 Array<DP,1> time_para,
 Array<DP,1> time_save,
 Array<int,1> misc_output_para,
 Array<int,1> ET_parameters,
 Array<DP,1> ET_shell_radii_sector_array,
 Array<int,1> no_output_k_r,
 Array<int,2> out_k_r_array,
 Array<int,1> init_cond_meta_para,
 Array<int,1> init_cond_int_para,
 Array<DP,1> init_cond_double_para,
 string init_cond_string_para[],
 Array<int,1> force_meta_para,
 Array<int,1> force_int_para,
 Array<DP,1> force_double_para,
 string force_string_para[]
 )
{

	string s;
	int i;

	if (my_id == master_id)
		cout << endl << "*********** Field parameters ***********" << endl;

	// Every processor checks if the para_file can be opened.
	if (! para_file.is_open()) {
		cout << "Unable to open para_file; Exiting Program: MY_ID " << my_id << endl;
		exit(1);
	}

	//
	//	N[i]
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	if (my_id == master_id)
		cout << "N[i]: " ;
	N[0] = 0;

	for (i = 1; i<= dim; i++) {
		para_file >> N[i];
		if (my_id == master_id)
			cout << N[i] << "   " ;
	}
	if (my_id == master_id)
		cout << endl;

	//
	//	Solver-para
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	para_file >> s >> solver_meta_para(1) >> solver_meta_para(2) >> solver_meta_para(3);
	if (my_id == master_id)	{
		cout << endl << "*********** Solver parameters ***********" << endl;
		cout << "Solver meta para (intgers): "			<< solver_meta_para(1) << endl;
		cout << "Solver meta para (double): "			<< solver_meta_para(2) << endl;
		cout << "Solver meta para (string): "			<< solver_meta_para(3) << endl << endl;
	}
	if (my_id == master_id)
		cout << "Solver integer parameters: ";

	solver_int_para = 0;
	solver_double_para = 0.0;

	para_file >> s;
	for (i = 1; i <= solver_meta_para(1); i++) {
		para_file >> solver_int_para(i);
		if (my_id == master_id)
			cout << solver_int_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	if (my_id == master_id)
		cout << "Solver double parameters: ";

	para_file >> s;
	for (i = 1; i <= solver_meta_para(2); i++) {
		para_file >> solver_double_para(i);
		if (my_id == master_id)
			cout << solver_double_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	getline(para_file, s);

	if (my_id == master_id)
		cout << "Solver string parameters: ";

	para_file >> s;
	for (i = 1; i<= solver_meta_para(3); i++) {
		para_file >> solver_string_para[i];
		if (my_id == master_id)
			cout <<  solver_string_para[i];
	}
	if (my_id == master_id)
		cout << endl;

	getline(para_file, s); getline(para_file, s); getline(para_file, s);


	//
	// String_switches
	//
	for (i = 1; i <= 5 ; i++) {
		para_file >> s;
		para_file >> string_switches[i];
	}

	if (my_id == master_id)	{
		cout << endl << "*********** String switches **************" << endl;
		cout << "Basis_type: (FOUR or SCFT)"			<<  string_switches[1] << endl;
		cout << "Alias_switch: "				<<  string_switches[2] << endl;
		cout << "Integration scheme: "				<<  string_switches[3] << endl;
		cout << "number input field mode (ASCII or BINARY): "	<<  string_switches[4] << endl;
		cout << "number output field mode (ASCII or BINARY): "	<<  string_switches[5] << endl;
	}

	//
	// Switches
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	for (i = 1; i <= 22; i++) {
		para_file >> s;
		para_file >> switches(i);
	}
	if (my_id == master_id)	{
		cout << endl << "*********** Integer switches **************" << endl;
		cout << "Hyperdissipation switch: " 		<< switches(1) << endl;
		cout << "Low-dimension switch: "		<< switches(2) << endl;
		cout << "LES switch: "				<< switches(3) << endl;
		cout << "Anisotropic ring switch: "		<< switches(4) << endl;
		cout << "Anisotropic cylinder switch: "		<< switches(5) << endl;
		cout << "structure_fn_switch: "			<< switches(6) << endl;
		cout << "planar_structure_fn_switch: "		<< switches(7) << endl;
		cout << "skpq_switch: "				<< switches(8) << endl;
		cout << "output_real_space_switch: "		<< switches(9) << endl;
		cout << "free-slip-verticalwall-switch: "	<< switches(11) << endl;
		cout << "helicity-flux_switch: "		<< switches(12) << endl;
		cout << "helicity-shells_switch: "		<< switches(13) << endl;
		cout << "anisotropic-switch (1,2,3): "		<< switches(14) << endl;
		cout << "waveno-switch (0: actual, 1: grid): "  << switches(15) << endl;
		cout << "fixed_dt_switch: "			<< switches(16) << endl;
		cout << "apply_realitycond_IC_switch: "		<< switches(17) << endl;
		cout << "apply_realitycond_alltime_switch: "	<< switches(18) << endl;
		cout << "output_vx_vy_switch: "			<< switches(19) << endl;
		cout << "input_vx_vy_switch: "			<< switches(20) << endl;
		cout << "force_W_switch: "			<< switches(21) << endl;
		cout << "force_T_switch: "			<< switches(22) << endl;
	}

	//
	// Dissipation coefficients
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s); getline(para_file, s);
	if (my_id == master_id)
		cout << endl << "*********** Dissipation coefficients **************" << endl;
		cout << "Dissipation coefficients: " ;

	for (i = 0; i < number_of_fields; i++) {
		para_file >> diss_coefficient[i];
		if (my_id == master_id)
			cout << diss_coefficient[i] << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	//
	// Hyper-dissipation coefficients.  Default value = 0
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	if (my_id == master_id)
		cout << "Hyper dissipation coefficients: " ;

	for (i = 0; i < number_of_fields; i++) {
		para_file >> hyper_diss_coefficient[i];
		if (my_id == master_id)
			cout << hyper_diss_coefficient[i] << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	//
	//	Time-para and Time-save
	//
	getline(para_file, s); getline(para_file, s);  getline(para_file, s);
	if (my_id == master_id)
		cout << endl << "*********** Time parameters **************" << endl;
		cout << "time_para (t0, tfinal, dt_min): ";

	for (i = 1; i <= 5; i++) {
		para_file >> s;
		para_file >> time_para(i);

		if (my_id == master_id)
			cout << time_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	if (my_id == master_id)
		cout << "time_save: ";
	for (i = 1; i <= 19; i++) {
		para_file >> s;
		para_file >> time_save(i);
		if (my_id == master_id)
			cout << time_save(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	//
	// spectrum fn, structure fn, N_out_reduced etc.
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	for (i = 1; i <= 5; i++) {
		para_file >> s;
		para_file >> misc_output_para(i);
	}

	if (my_id == master_id)	{
		cout << endl << "*********** Spectrum parameters **************" << endl;
		cout << "Spectrum parameter: sector-input-scheme: "<< misc_output_para(1) << endl;
		cout << "Spectrum parameter: no_sector-spectrum: " << misc_output_para(2) << endl;
		cout << "Spectrum parameter: no_cylinder_slabs_spectrum " << misc_output_para(3) << endl;
		cout << "Structure fn parameter: q-min " << misc_output_para(4) << endl;
		cout << "Structure fn parameter: q-max " << misc_output_para(5) << endl;
	}

	//
	//  N_out_reduced
	//
	if (my_id == master_id)
		cout << "N_out_reduced[i]: " ;

	para_file >> s;
	for (i = 6; i<= 5+dim; i++) {
		para_file >> misc_output_para(i);
		cout << misc_output_para(i) << "   " ;
	}
	if (my_id == master_id)
		cout << endl;

	//
	// Energy transfer parameters
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	for (i = 1; i < 12; i++) {
		para_file >> s;
		para_file >> ET_parameters(i);
	}

	if (my_id == master_id) {
		cout << endl << "*********** Energy transfer parameters **************" << endl;
		cout << "Energy transfer (real_imag_switch, nospheres, noshells):  "
		     << ET_parameters(1) << "   "  << ET_parameters(2) << "   "
		     << ET_parameters(3) << endl;
		cout << "Rings (no-ring-shells, no-sectors): "
		     << ET_parameters(4) << "   "  << ET_parameters(5)  << endl;
		cout << "Cylinders (no-cylinder-shells, no-cylinder-kpll-planes): "
		     << ET_parameters(6) << "   "  << ET_parameters(7)  << endl;
		cout << "Input-schemes (shell, ring-shell, sector, cylinder-shell, cylinder-kpll): "
		     << ET_parameters(8) << "   "  << ET_parameters(9)  << "  "
		     << ET_parameters(10) << "   "  << ET_parameters(11)  << "  "
		     << ET_parameters(12) << endl;
	}

	// At present, reading of shells, sectors etc. not implemented.
	ET_shell_radii_sector_array = 0.0;


	//
	// read coordinates of wavenumber and positions at which the field are to be printed.
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	if (my_id == master_id)
		cout << endl << "*********** Wavenumber positions **************" << endl;
	para_file >> s >> no_output_k_r(1);
	if (my_id == master_id)
		cout << "No of fourier-space points for which field values and Tk are to be printed: "
		     << no_output_k_r(1) << endl;

	para_file >> s >> no_output_k_r(2);

	if (my_id == master_id)
		cout << "No of real-space points for which the field are to be printed: "
		     << no_output_k_r(2) << endl;


	// read \vec{k}
	//
	getline(para_file, s); getline(para_file, s);
	int kx, ky, kz;			// wavenumbers
	out_k_r_array = 0;    // initialize to zero.

	for (i=1; i <= no_output_k_r(1); i++) {
		para_file >> kx >> ky >> kz;

		out_k_r_array(i,1) = kx;
		out_k_r_array(i,2) = ky;
		out_k_r_array(i,3) = kz;

		if (my_id == master_id)
			cout << "kx, ky, kz: " <<  kx <<"   " << ky << "   "  << kz << endl;
	}


	// read \vec{r}
	//
	getline(para_file, s); getline(para_file, s);
	for (i=no_output_k_r(1)+1; i<=no_output_k_r(1)+no_output_k_r(2); i++) {
		para_file >> kx >> ky >> kz;

		out_k_r_array(i,1) = kx;
		out_k_r_array(i,2) = ky;
		out_k_r_array(i,3) = kz;

		if (my_id == master_id)
			cout << "rx, ry, rz: " <<  kx <<"   " << ky << "   "  << kz << endl;
	}


	//
	//	Initial condition
	//
	getline(para_file, s); getline(para_file, s); getline(para_file, s);
	if (my_id == master_id)
		cout << endl << "*********** Initial condition parameters ***********" << endl;

	para_file >> s >> init_cond_meta_para(1) >> init_cond_meta_para(2) >> init_cond_meta_para(3);
	if (my_id == master_id) {
		cout << "Init cond meta para (intgers): " << init_cond_meta_para(1) << endl;
		cout << "Init cond meta para (double): "  << init_cond_meta_para(2) << endl;
		cout << "Init cond meta para (string): "  << init_cond_meta_para(3) << endl;
	}

	init_cond_int_para = 0;
	init_cond_double_para = 0.0;

	if (my_id == master_id)
		cout << "Init cond integer parameters: ";
	para_file >> s;
	for (i = 1; i <= init_cond_meta_para(1); i++) {
		para_file >> init_cond_int_para(i);
		if (my_id == master_id)
			cout << init_cond_int_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	if (my_id == master_id)
		cout << "Init cond double parameters: ";
	para_file >> s;
	for (i = 1; i <= init_cond_meta_para(2); i++) {
		para_file >> init_cond_double_para(i);
		if (my_id == master_id)	cout << init_cond_double_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;
	getline(para_file, s);

	if (my_id == master_id)
		cout << "Init cond string parameters: ";
	para_file >> s;
	for (i=1; i<=init_cond_meta_para(3); i++) {
		para_file >> init_cond_string_para[i];
		if (my_id == master_id)		cout <<  init_cond_string_para[i] << endl;
	}

	getline(para_file, s); getline(para_file, s); getline(para_file, s);

	//
	//	Forcing
	//
	if (my_id == master_id)
		cout << endl << "*********** Forcing parameters ***********" << endl;

	para_file >> s >> force_meta_para(1) >> force_meta_para(2) >> force_meta_para(3);

	if (my_id == master_id) {
		cout << "Force meta para (intgers): " << force_meta_para(1) << endl;
		cout << "Force meta para (double): "  << force_meta_para(2) << endl;
		cout << "Force meta para (string): "  << force_meta_para(3) << endl;
	}

	force_int_para = 0;
	force_double_para = 0.0;

	if (my_id == master_id)
		cout << "Force integer parameters: ";
	para_file >> s;
	for (i = 1; i <= force_meta_para(1); i++) {
		para_file >> force_int_para(i);
		if (my_id == master_id)
			cout << force_int_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	if (my_id == master_id)
		cout << "Force double parameters: ";
	para_file >> s;
	for (i = 1; i <= force_meta_para(2); i++) {
		para_file >> force_double_para(i);
		if (my_id == master_id)
			cout << force_double_para(i) << "    ";
	}
	if (my_id == master_id)
		cout << endl;

	getline(para_file, s);

	if (my_id == master_id)
		cout << "Force string parameters: ";
	para_file >> s;
	for (i=1; i<=force_meta_para(3); i++) {
		para_file >> force_string_para[i];
		if (my_id == master_id)	cout << force_string_para[i];
	}
	if (my_id == master_id)
		cout << endl;

	getline(para_file, s); getline(para_file, s); getline(para_file, s);
}
