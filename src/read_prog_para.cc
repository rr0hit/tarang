/*! \file  read_prog_para.cc
 *
 *	@brief  Reads prog prameters: prog_para_file, prog_kind, data_dir_name, D.
 *
 *	@author  M. K. Verma
 *	@version 4.0  MPI
 *	@date Sept 2008
 */


#include "main.h"

void Read_prog_para(string &kind, string& data_dir_name)
{
	ifstream prog_para_file;				// prog_para_file defined only in the master node
	string str;

	prog_para_file.open("prog_para.d");

	if (my_id == master_id)
		cout << endl << "======= Reading program parameters =========" << endl << endl << endl;

	if (! prog_para_file.is_open())
		{
			cout << "Unable to open prog_para_file: MY_ID  " << my_id <<  endl;
			exit(1);
		}

	while (!getline(prog_para_file, str).eof()) {
		if (!str.length() || str[0] == '#')
			continue; /* Ignore comments, blank lines */
		else if (str[0] == '$') {
			str.erase(0, 1); /* Erase the leading '$' */
			kind = str;
		} else
			data_dir_name = str;
	}

	if (my_id == master_id) {
		cout << "Program: " << kind << endl;
		cout << "data_dir_name: " << data_dir_name << endl;
		cout  << "===========================================" << endl << endl;
	}
}
