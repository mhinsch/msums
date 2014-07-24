#ifndef MSDATAFILE_H
#define MSDATAFILE_H


#include <iostream>
#include <string>
#include <vector>

#include "sample.h"

using namespace std;


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
void read_dataset(istream & inp, 
	size_t dataset, 
	const vector<vector<size_t> > & nSeq, const vector<size_t> & nSites, 
	vector<vector<Sample<string, char> > > & seqhs, vector<string> & seqOls);


#endif	// MSDATAFILE_H
