#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "spmain.h"
#include "sputil.h"
#include "spinout.h"

using namespace std;
using namespace boost;

typedef SPIOException SPIOE;

void skip_space(istream & inp, string & str)
	{
	while(getline(inp, str) && all(str, is_space()));
	}

void get_initial_conditions
	(const string & inputfile, vector<vector<int> > & nSeq, vector<int> & nSites, 
	 int & ndatasets, string & datafile)
 	{
/* reads inputfile for values of parameters and check them against legal range*/
	ifstream inp(inputfile.c_str());
	if (!inp.good())
		throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot open file " + 
			inputfile);
	
	int nloci, npops;

	if (!get_value(inp, nloci))
		throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read nloci");
	if (!get_value(inp, npops))
	   throw SPIOE(ERR_LOC
			" Error in reading initial conditions: cannot read npops");

	//nSeq.reserve(npops);

	for(int p=0; p<npops; p++)
		{
		nSeq.push_back(vector<int>());
		nSeq.back().resize(nloci);
		}

	nSites.resize(nloci);

	for (int l=0; l<nloci; l++)
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
			nSeq[p][l] = v;
			}

		if (!get_value(inp, nSites[l]))
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
	@param nSeq #sequences[population][locus].
	@param nSites number of sites[locus].
	@param seqhs sequences[population][locus][number].
	@param seqOls outgroup data[locus].
	@param nspolyl 
	*/
void get_dataset(istream & inp, 
	long dataset, const vector<vector<int> > & nSeq, const vector<int> & nSites, 
	vector<vector<vector<string> > > & seqhs, vector<string> & seqOls, 
	vector<int> & nspolyl)
	{
	int nseg;

	const string match_segsites = "segsites:";
	const string match_positions = "positions:";
	string str;

	const int nPops = nSeq.size();
	const int nloci = nSites.size();

	// loop over all loci
	for (int l=0; l<nloci; l++)
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

		nspolyl[l] = nseg;

		if (nseg > nSites[l]) 
			{	/*if too many segregating sites*/
			nspolyl[l] = nSites[l];
			// TODO: ask Ludovic whether error message in this case
			//errorfile << "\nerror in reading dataset: nseg=" << nseg 
			//	<< " > total number of sites=" << nSites[l] 
			//	<< " in locus " << l;
			}	/*of else*/

		skip_space(inp, str); // skip blanks and 'positions:' line

		// for all species
		for (int p=0; p<nPops; p++)
			{
			// all sequences of species p
			for (int h=0; h<nSeq[p][l]; h++) 
				{	
				// fill with nSites '0's
				//seqhs[p][l][h].insert(0, nSites[l], '0');
				seqhs[p][l][h].clear();
				seqhs[p][l][h].resize(nSites[l], '0');

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

				seqhs[p][l][h].replace(0, str.size(), str);
				}
			}	/* end loop over haplotypes in species A*/
			
		seqOls[l].clear();
		seqOls[l].resize(nSites[l], '0');
		//seqOls[l].insert(0, nSites[l], '0');
		}	/*end loop over loci*/
	}		 /*end of get_dataset*/


/*****************************************************************************/

void print_initial_conditions(ostream & out, 
	const string & inputfile, vector<vector<int> > nSeq, vector<int> nSites, 
	long ndatasets, const string & datafile) 
/* print initial values of parameters read from input file as a check*/
	{
	const int nLoci = nSites.size();

	out << "\nThe following informations have been successfully loaded from " 
		<< inputfile << ":";
	out << "\n\tnumber of populations: " << nSeq.size();
	out << "\n\tnumber of loci: " << nLoci;
	for(int l=0;l<nLoci;l++)	/*loop over loci*/
		{
		for (size_t p=0; p<nSeq.size(); p++)
			out << "\n\tnumber of sequences in species " << p << 
				" at locus " << l << ":" << nSeq[p][l];

		out << "\n\t\tnumber of sites at locus " << l << ": " << nSites[l];
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



