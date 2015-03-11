#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "msdatafile.h"
#include "sputil.h"

using namespace std;

typedef SPIOException SPIOE;


/** Reads datasetfile and extracts DNA sequences from all population and 
	the outgroup.
	@param inp file to read data from.
	@param dataset number of current data set.
	@param n_sequences #sequences[population][locus].
	@param n_sites number of sites[locus].
	@param sequences sequences[population][locus][number].
	@param outgroup outgroup data[locus].
	@param mask missing data bool[locus][site].
	*/
void read_dataset(istream & inp, size_t dataset, 
	const vector<vector<size_t> > & n_sequences, const vector<size_t> & n_sites, 
	vector<vector<StrSample> > & sequences, vector<Sequence> & outgroups,
	const vector<vector<bool> > & mask )
	{
	size_t nseg;

	const string match_segsites = "segsites:";

	vector<int> positions;

	string str;

	// this is a bit non-obvious
	const size_t n_pops = n_sequences.size();
	const size_t n_loci = n_sites.size();

	// loop over all loci
	for (size_t l=0; l<n_loci; l++)
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

		if (nseg > n_sites[l]) 
			{	/*if too many segregating sites*/
			nseg = n_sites[l];
			// TODO: ask Ludovic whether error message in this case
			//errorfile << "\nerror in reading dataset: nseg=" << nseg 
			//	<< " > total number of sites=" << n_sites[l] 
			//	<< " in locus " << l;
			}

		skip_space(inp, str); // next non-blank line
		if (mask.size())	// we need position info
			{
			positions.clear();
			positions.resize(nseg);
			
			itmp.clear();
			itmp.str(str);
			// throw away first part of line
			stmp = "";
			itmp >> stmp;

			float p;
			for (size_t s=0; s<nseg; s++)
				{
				itmp >> p;
				positions[s] = lround(p * mask[l].size());
				}
			}

		// for all population
		for (size_t p=0; p<n_pops; p++)
			{
			sequences[p][l].set_tot_n_sites(n_sites[l]);

			// all sequences of population p
			for (size_t h=0; h<n_sequences[p][l]; h++) 
				{	
				string str("");
				str.reserve(nseg);

				// careful, this is two things in one
				if (nseg>0 && !getline(inp, str))
					throw SPIOE(
						string(ERR_LOC 
							" error in reading dataset: cannot read seq of "
							"haplotype ") + 
						lexical_cast<string>(h) + " in population " +
						lexical_cast<string>(p) + " at locus " +
						lexical_cast<string>(l));

				Sequence & seq = sequences[p][l].sequence(h);
				seq.clear(); seq.reserve(nseg);
				if (mask.size())
					{
					for (size_t i=0; i<str.size(); i++)
						// only add non-masked sites
						if (mask[l][positions[i]] == 0)
							seq.push_back(int(str[i]) - int('0'));
					}
				else
					for (size_t i=0; i<str.size(); i++)
						{
						if (str[i] != '0' && str[i] != '1')
							cerr << "unexpected input: " << str[i] << std::endl;
						seq.push_back(int(str[i]) - int('0'));
						}
				}
			}	/* end loop over haplotypes in population A*/
			
		outgroups[l].clear();
		outgroups[l].resize(nseg, 0);
		}	/*end loop over loci*/
	}		 /*end of get_dataset*/


