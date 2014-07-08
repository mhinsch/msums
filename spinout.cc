#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "spmain.h"
#include "sputil.h"
#include "spinout.h"

using namespace std;

typedef SPIOException SPIOE;

void skip_space(istream & inp, string & str)
	{
	while(getline(inp, str) && boost::all(str, boost::is_space()));
	}

void get_initial_conditions
	(const string & inputfile, 
	 vector<vector<int> > & n_sequences, vector<int> & n_sites, 
	 int & ndatasets, string & datafile)
 	{
/* reads inputfile for values of parameters and check them against legal range*/
	ifstream inp(inputfile.c_str());
	if (!inp.good())
		throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot open file " + 
			inputfile);
	
	int n_loci, npops;

	if (!get_value(inp, n_loci))
		throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read nloci");
	if (!get_value(inp, npops))
	   throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read npops");

	//n_sequences.reserve(npops);

	for(int p=0; p<npops; p++)
		{
		n_sequences.push_back(vector<int>());
		n_sequences.back().resize(n_loci);
		}

	n_sites.resize(n_loci);

	for (int l=0; l<n_loci; l++)
		{
		for (int p=0; p<npops; p++)
			{
			int v;
			if (!get_value<int>(inp, v))
				{
				throw SPIOE(ERR_LOC
					" Error in reading initial conditions: cannot read nSeq at locus "
						+ lexical_cast<string>(l));
				}
			n_sequences[p][l] = v;
			}

		if (!get_value(inp, n_sites[l]))
			{
			throw SPIOE(ERR_LOC
				" Error in reading initial conditions: cannot read nSites at locus "
					+ lexical_cast<string>(l));
			}
		}

	if (!get_value(inp, ndatasets))
		throw SPIOE(ERR_LOC
			" Error in reading inital conditions: cannot read ndatasets");
	if (!get_value(inp, datafile))
	   throw SPIOE(ERR_LOC	
			" Error in reading initial conditions: cannot read data file name");
	} /*end of get_initial_conditions*/

/*****************************************************************************/

/** Reads datasetfile and extracts DNA sequences from all species and 
	the outgroup.
	@param inp file to read data from.
	@param dataset number of current data set.
	@param n_sequences #sequences[population][locus].
	@param n_sites number of sites[locus].
	@param sequences sequences[population][locus][number].
	@param outgroup outgroup data[locus].
	@param n_polym_sites 
	*/
void get_dataset(istream & inp, 
	long dataset, 
	const vector<vector<int> > & n_sequences, const vector<int> & n_sites, 
	vector<vector<Sample> > & sequences, Sample & outgroup, 
	vector<int> & n_polym_sites)
	{
	int nseg;

	const string match_segsites = "segsites:";
	const string match_positions = "positions:";
	string str;

	// this is a bit non-obvious
	const int n_pops = n_sequences.size();
	const int n_loci = n_sites.size();

	// loop over all loci
	for (int l=0; l<n_loci; l++)
		{
		istringstream itmp;
		string stmp;
		// skip all lines that don't start with 'segsites:'
		do
			{
			// next non-empty line
			skip_space(inp, str);
			itmp.clear();
			itmp.str(str);
			// get first part of line
			stmp = "";
			itmp >> stmp;

			if (!inp.good() || itmp.fail())
				throw SPIOE(string(ERR_LOC " Error in reading dataset ") +
					lexical_cast<string>(dataset) + " at locus " +
					lexical_cast<string>(l) + ": expected '" + match_segsites +
					"', found '" + stmp + "'");
			}
		while (stmp != match_segsites);

		itmp >> nseg;

		if (itmp.fail())
			throw SPIOE(string(ERR_LOC " Error in reading dataset ") + 
				lexical_cast<string>(dataset) + " at locus " + 
				lexical_cast<string>(l) +  ": couldn't read number of seg sites ");

		n_polym_sites[l] = nseg;

		if (nseg > n_sites[l]) 
			{	/*if too many segregating sites*/
			n_polym_sites[l] = n_sites[l];
			// TODO: ask Ludovic whether error message in this case
			//errorfile << "\nerror in reading dataset: nseg=" << nseg 
			//	<< " > total number of sites=" << n_sites[l] 
			//	<< " in locus " << l;
			}

		skip_space(inp, str); // skip blanks and 'positions:' line

		// for all species
		for (int p=0; p<n_pops; p++)
			{
			// all sequences of species p
			for (int h=0; h<n_sequences[p][l]; h++) 
				{	
				// fill with n_sites '0's
				//sequences[p][l][h].insert(0, n_sites[l], '0');
				sequences[p][l][h].clear();
				sequences[p][l][h].resize(n_sites[l], '0');

				if (nseg == 0)	// no sites, skip reading
					continue;

				if (!getline(inp, str))
					throw SPIOE(
						string(ERR_LOC 
							" error in reading dataset: cannot read seq of "
							"haplotype ") + 
						lexical_cast<string>(h) + " in species " +
						lexical_cast<string>(p) + " at locus " +
						lexical_cast<string>(l));

				sequences[p][l][h].replace(0, str.size(), str);
				}
			}	/* end loop over haplotypes in species A*/
			
		outgroup[l].clear();
		outgroup[l].resize(n_sites[l], '0');
		//outgroup[l].insert(0, n_sites[l], '0');
		}	/*end loop over loci*/
	}		 /*end of get_dataset*/



/*****************************************************************************/

/* print initial values of parameters read from input file as a check*/
void print_initial_conditions(ostream & out, 
	const string & inputfile, 
	vector<vector<int> > n_sequences, vector<int> n_sites, 
	long ndatasets, const string & datafile) 
	{
	const int nLoci = n_sites.size();

	out << "\nThe following informations have been successfully loaded from " 
		<< inputfile << ":";
	out << "\n\tnumber of populations: " << n_sequences.size();
	out << "\n\tnumber of loci: " << nLoci;
	for(int l=0;l<nLoci;l++)	/*loop over loci*/
		{
		for (size_t p=0; p<n_sequences.size(); p++)
			out << "\n\tnumber of sequences in species " << p << 
				" at locus " << l << ":" << n_sequences[p][l];

		out << "\n\t\tnumber of sites at locus " << l << ": " << n_sites[l];
		}

	out << "\n\tnumber of replicate datasets: " << ndatasets; 
	out << "\n\tgeneric name of the datasets: " << datafile;
	out << "\n\n";
	}	  /*end of print_initial_conditions*/


/*****************************************************************************/

/** Initialize the output file where writing statistics of shared and fixed 
	polymorphism to outputfile one line for each replicate dataset. */
void initialize_write_ABCstat(ostream & file, 
	const vector<bool> & pops, const vector<vector<bool> > & pairs, 
	const vector<bool> sstats, const vector<bool> pstats,
	const StatDict & dict, bool writeAggr)
	{
	// TODO: labels according to ludovic's specs
	//
	file << "\ndataset";
	
	// works for now, might have to be separate param in the future
	if (!writeAggr)
		file << "\tlocus";

	for (size_t p=0; p<pops.size(); p++)
		{
		if (!pops[p])
			continue;
		
		for (size_t s=0; s<sstats.size(); s++)
			{
			// not shown
			if (!sstats[s])
				continue;

			file << '\t' << dict.singles[s] << '_' << p;
			if (writeAggr)
				{
				file << "_mean";
				file << '\t' << dict.singles[s] << '_' << p << "_std";
				}
			}
		}

	for (size_t p=0; p<pairs.size()-1; p++)
		for (size_t pp=p+1; pp<pairs.size(); pp++)
			{
			if (!pairs[p][pp])
				continue;

			for (size_t s=0; s<pstats.size(); s++)
				{
				if (!pstats[s])
					continue;
				
				file << '\t' << dict.pairs[s] << '_' << p << 'x' << pp;
				if (writeAggr)
					{
					file << "_mean";
					file << '\t' << dict.pairs[s] << '_' << p << 'x' << 
						pp << "_std";
					}
				}
			}

	file << endl;
	}	/*end of initialize_write_ABCstat*/



/*****************************************************************************/

/*write statistics of shared and fixed polymorphism to outputfile one line for each replicate dataset*/
//void write_ABCstat(ostream & out, 
//	long dataset, int nloc, const All_stats &r)
//	{
//	out << "\n" << dataset;
//	out << "\t" << r.pearson_corr_pi();
//	}	/*end of write_ABCstat*/



