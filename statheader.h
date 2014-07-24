#ifndef STATHEADER_H
#define STATHEADER_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;


/** Initialize the output file where writing statistics of shared and fixed 
	polymorphism to outputfile one line for each replicate dataset. */
void write_ABCstat_header(ostream & file, 
	const vector<bool> & pops, const vector<vector<bool> > & pairs, 
	const vector<string> & ss_names, const vector<string> & ps_names,
	bool writeAggr);



#endif	//STATHEADER_H
