#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "sputil.h"
#include "configfile.h"

using namespace std;

typedef SPIOException SPIOE;

void ConfigFile::read(istream & inp)
 	{
	int n_loci, npops;

	if (!get_value(inp, n_loci))
		throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read nloci");
	if (!get_value(inp, npops))
	   throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read npops");

	for(size_t p=0; p<npops; p++)
		{
		_n_sequences.push_back(vector<size_t>());
		_n_sequences.back().resize(n_loci);
		}

	_n_sites.resize(n_loci);

	for (size_t l=0; l<n_loci; l++)
		{
		for (size_t  p=0; p<npops; p++)
			{
			int v;
			if (!get_value<int>(inp, v))
				{
				throw SPIOE(ERR_LOC
					" Error in reading initial conditions: cannot read nSeq at locus "
						+ lexical_cast<string>(l));
				}
			_n_sequences[p][l] = v;
			}

		if (!get_value(inp, _n_sites[l]))
			{
			throw SPIOE(ERR_LOC
				" Error in reading initial conditions: cannot read nSites at locus "
					+ lexical_cast<string>(l));
			}
		}

	if (!get_value(inp, _ndatasets))
		throw SPIOE(ERR_LOC
			" Error in reading inital conditions: cannot read ndatasets");
	if (!get_value(inp, _datafile))
	   throw SPIOE(ERR_LOC	
			" Error in reading initial conditions: cannot read data file name");
	} 


/* print initial values of parameters read from input file as a check*/
void ConfigFile::dump(ostream & out)
	{
	const size_t nLoci = _n_sites.size();

	out << "\n\tnumber of populations: " << _n_sequences.size();
	out << "\n\tnumber of loci: " << nLoci;
	for(int l=0;l<nLoci;l++)	/*loop over loci*/
		{
		for (size_t p=0; p<_n_sequences.size(); p++)
			out << "\n\tnumber of sequences in population " << p << 
				" at locus " << l << ":" << _n_sequences[p][l];

		out << "\n\t\tnumber of sites at locus " << l << ": " << _n_sites[l];
		}

	out << "\n\tnumber of replicate datasets: " << _ndatasets; 
	out << "\n\tname of the dataset file: " << _datafile;
	out << "\n\n";
	}	  /*end of print_initial_conditions*/

