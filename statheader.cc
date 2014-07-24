#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "statheader.h"


using namespace std;

void write_stat_name(ostream & file, const string & name, 
	bool aggr = false, const string & sep = "\t")
	{
	file << sep << name;
	if (aggr)
		{
		file << "_mean";
		file << sep << name << "_std";
		}
	}

string stat_name(const string & name, size_t p)
	{
	ostringstream oss(name);
	oss << '_' << p;
	return oss.str();
	}
	
string stat_name(const string & name, size_t p1, size_t p2)
	{
	ostringstream oss(name);
	oss << '_' << p1 << 'x' << p2;
	return oss.str();
	}

string stat_name(const string & name, const vector<size_t> & group)
	{
	ostringstream oss(name);
	oss << '_';
	copy(group.begin(), group.end(), ostream_iterator<size_t>(oss, "x"));
	return oss.str();
	}

/*****************************************************************************/

/** Initialize the output file where writing statistics of shared and fixed 
	polymorphism to outputfile one line for each replicate dataset. */
void write_ABCstat_header(ostream & file, 
	const vector<bool> & pops, const vector<vector<bool> > & pairs, 
	const vector<string> & ss_names, const vector<string> & ps_names,
	const vector<string> & ms_names, const vector<vector<size_t> > & groups,
	bool write_aggr)
	{
	file << "\ndataset";
	
	// works for now, might have to be separate param in the future
	if (!write_aggr)
		file << "\tlocus";

	for (size_t p=0; p<pops.size(); p++)
		{
		if (!pops[p])
			continue;
		
		for (size_t s=0; s<ss_names.size(); s++)
			write_stat_name(file, stat_name(ss_names[s], p), write_aggr);
		}

	for (size_t p=0; p<pairs.size()-1; p++)
		for (size_t pp=p+1; pp<pairs.size(); pp++)
			{
			if (!pairs[p][pp])
				continue;

			for (size_t s=0; s<ps_names.size(); s++)
				write_stat_name(file, stat_name(ps_names[s], p, pp), write_aggr);
			}

	for (size_t s=0; s<ms_names.size(); s++)
		write_stat_name(file, stat_name(ms_names[s], groups[s]), write_aggr);

	file << endl;
	}



