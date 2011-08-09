/*! \file  main.h
 *
 *	@brief  Main file for Tarang (MPI). Calls various modules.
 *
 *	@author  M. K. Verma
 *	@version 4.0  MPI
 *	@date Sept 2008
 */

#include "IncFluid.h"
// Contains IncFlow declarations (which contains fields +fourier+sincosfour defs)

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>

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
	Array<int,1> input_meta_para,
	Array<int,1> input_int_para,
	Array<DP,1> input_double_para,
	string input_string_para[],
	Array<int,1> force_meta_para,
	Array<int,1> force_int_para,
	Array<DP,1> force_double_para,
	string force_string_para[]
	);

void Read_diag_para
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
	Array<int,1>  diagnostic_procedure,
	Array<DP,1> time_para,
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
	);

void Read_prog_para(string &kind, string& data_dir_name);

int Ifluid_main(string data_dir_name);
int Ifluid_diag_main(string data_dir_name);
int Iscalar_main(string data_dir_name);
int Iscalar_diag_main(string data_dir_name);
int IMHD_main(string data_dir_name);
int IMHD_diag_main(string data_dir_name);
int RB_slip_main(string data_dir_name);
int RB_slip_diag_main(string data_dir_name);
int RB_slipMHD_main(string data_dir_name);
int NonBoussinesq_main(string data_dir_name);
