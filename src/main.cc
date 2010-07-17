/* Tarang-4.0
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-4.0 .
 *
 * Tarang-4.0 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-4.0 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-4.0; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */


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

fftw_plan r2c_plan_FOUR, c2r_plan_FOUR;
fftw_plan r2c_plan_SCFT, c2r_plan_SCFT, sintr_plan_SCFT, costr_plan_SCFT, 
					isintr_plan_SCFT, icostr_plan_SCFT;  // i for inverse

fftw_plan r2c_1d_plan_SCFT, c2r_1d_plan_SCFT;  // for 2D
					
Uniform<DP> SPECrand;					// Global variable for random no generation		
			
					
//********************************************************************************************* 


int main(int argc, char** argv)
{

	string		prog_kind;
	string		data_dir_name;					// data is stored here
	
	time_t		start;
	time_t		end;
	
	// start main program
	
	SPECrand.seed((unsigned int)time(0));				// Initialize random number seed
	
	// MPI INIT
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	time (&start);
	
	// All procs read program parameters

	ifstream prog_para_file;				// prog_para_file defined only in the master node
	prog_para_file.open("prog_para.d");
	Read_prog_para(prog_para_file, prog_kind, data_dir_name); 
		
	if (prog_kind == "INC_FLUID")
		Ifluid_main(data_dir_name);

	else if (prog_kind == "INC_PASSIVE_SCALAR")
		Ipassive_scalar_main(data_dir_name);
		
	else if (prog_kind == "INC_MHD")
		IMHD_main(data_dir_name);
				
	else if (prog_kind == "RB_SLIP")
		RB_slip_main(data_dir_name);
		
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



//********************************** End of main.cc *******************************************


