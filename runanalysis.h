#ifndef RUNANALYSIS_H
#define RUNANALYSIS_H

#include <iostream>
#include <vector>
#include <string>

#include "anafunctors.h"
#include "sample.h"
#include "pairsample.h"

using namespace std;

typedef AnalysisBase<Sample<string, char> > single_ana_t;
typedef AnalysisBase<PairSample<string, char> > pair_ana_t;
typedef AnalysisBase<vector<Sample<string, char>*> > group_ana_t;

void analyse(ostream & abcstatfile, int dataset,
	vector<vector<single_ana_t::sample_t> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<single_ana_t*> & single_stats,
	const vector<pair_ana_t*> & pair_stats,
	const vector<group_ana_t*> & group_stats,
	const vector<vector<size_t> > & groups);

void analyse_aggr(ostream & abcstatfile, int dataset,
	vector<vector<single_ana_t::sample_t> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<single_ana_t*> & single_stats,
	const vector<pair_ana_t*> & pair_stats,
	const vector<group_ana_t*> & group_stats,
	const vector<vector<size_t> > & groups);

#endif	// RUNANALYSIS_H
