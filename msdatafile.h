#ifndef MSDATAFILE_H
#define MSDATAFILE_H


#include <iostream>
#include <string>
#include <vector>

#include "sample.h"
#include "setup.h"

using namespace std;


/** Reads datasetfile and extracts DNA sequences from all population and 
	the outgroup.
	@param inp file to read data from.
	@param dataset number of current data set.
	@param n_sequences #sequences[population][locus].
	@param n_sites number of sites[locus].
	@param sequences sequences[population][locus][number].
	@param outgroup outgroup data[locus].
	@param mask missing data bool_sequence[locus].
	*/
void read_dataset(istream & inp, 
	size_t dataset, 
	const vector<vector<size_t> > & nSeq, const vector<size_t> & nSites, 
	vector<vector<StrSample> > & seqhs, vector<Sequence> & seqOls,
	const vector<vector<bool> > & mask = vector<vector<bool> > ());


#endif	// MSDATAFILE_H
