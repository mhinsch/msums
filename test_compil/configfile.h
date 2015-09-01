#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <vector>
#include <iostream>
#include <string>

using namespace std;

class ConfigFile
	{
protected:
	string _datafilename;
	/** Number of sequences per population x locus.*/
	vector<vector<size_t> > _n_sequences;
	/** Number of sites per locus. Assumes all sequences per locus are the
	    same length. */
	vector<size_t> _n_sites;
	size_t _n_datasets;

public:
	void read(istream & in);
	void dump(ostream & out);

	size_t n_pops() const
		{
		return _n_sequences.size();
		}
	size_t n_loci() const
		{
		return _n_sites.size();
		}

	size_t n_datasets() const
		{
		return _n_datasets;
		}

	const vector<vector<size_t> > & n_sequences() const
		{
		return _n_sequences;
		}

	const vector<size_t> & n_sites() const
		{
		return _n_sites;
		}

	const string & datafilename() const
		{
		return _datafilename;
		}
	};


#endif	// CONFIGFILE_H
