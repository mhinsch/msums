#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#include <vector>
#include <iostream>
#include <string>

class ConfigFile
	{
protected:
	std::string _datafilename;
	/** Number of sequences per population x locus.*/
	std::vector<std::vector<size_t> > _n_sequences;
	/** Number of sites per locus. Assumes all sequences per locus are the
	    same length. */
	std::vector<size_t> _n_sites;
	size_t _n_datasets;

public:
	void read(std::istream & in);
	void dump(std::ostream & out);

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

	const std::vector<std::vector<size_t> > & n_sequences() const
		{
		return _n_sequences;
		}

	const string & datafilename() const
		{
		return _datafilename;
		}
	};


#endif	// CONFIGFILE_H
